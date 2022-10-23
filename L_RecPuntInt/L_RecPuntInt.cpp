#include "L_RecPuntInt.h"

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
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Display.H>
#include <FL/fl_ask.H>
#include <FL/Fl_File_Chooser.H>
#endif // __COMPAT_FLTK__



void L_TransfSimil2DCelda::swap(L_TransfSimil2DCelda &other)
{
	L_PUSH_EXECUTING_FN("L_TransfSimil2DCelda::swap");
	other.L_TransfSimil2DCelda_POD::swap(*this);
	other.calL.swap(calL);
	L_POP_EXECUTING_FN("L_TransfSimil2DCelda::swap");
}

L_TransfSimil2DCeldaLista *L_TransfSimil2DCeldaHash::buscaCelda(int i, int j, int k, int z, int nobj, L_TransfSimil2DCeldaNodo **bin)
{
	L_PUSH_EXECUTING_FN("L_TransfSimil2DCeldaHash::buscaCelda");
	L_TransfSimil2DCeldaLista *binL;
	L_TransfSimil2DCeldaNodo *binPtr;
	int nh;
	nh=(i+j*11L+k*111L+z*1111L+nobj*11111L)%nTabla;
	if (nh<0)
		nh+=nTabla;
	binL=&tabla[nh];
	for (binPtr=binL->root; binPtr!=NULL; binPtr=binPtr->sig)
	{
		if (i==binPtr->c.i && j==binPtr->c.j && k==binPtr->c.k && z==binPtr->c.z && nobj==binPtr->c.nobj)
			break;
	}
	*bin=binPtr;
	L_POP_EXECUTING_FN("L_TransfSimil2DCeldaHash::buscaCelda");
	return binL;
}

void L_TransfSimil2DCeldaHash::minivoto(int i, int j, int k, int z, int nobj, L_Calce &cal, L_List<L_Calce> *cal_pool)
{
	L_PUSH_EXECUTING_FN("L_TransfSimil2DCeldaHash::minivoto");
	L_TransfSimil2DCeldaNodo *bin;
	L_TransfSimil2DCeldaLista *binL;
	binL=buscaCelda(i,j,k,z,nobj,&bin);
    if (bin!=NULL)
	{
		bin->c.calL.push_back(cal, cal_pool);
	}
	else
	{
		L_TransfSimil2DCelda binTemp;
		binTemp.i=i;
		binTemp.j=j;
		binTemp.k=k;
		binTemp.z=z;
		binTemp.nobj=nobj;
		binTemp.calL.push_back(cal, cal_pool);
		binL->push_back_swapping(binTemp);
	}
	L_POP_EXECUTING_FN("L_TransfSimil2DCeldaHash::minivoto");
}

void L_TransfSimil2DCeldaHash::swap(L_TransfSimil2DCeldaHash &other)
{
	L_TransfSimil2DCeldaLista *tabla_t;

	other.L_TransfSimil2DCeldaHash_POD::swap(*this);
	paramTSCH.swap(other.paramTSCH);
	//
	tabla_t = other.tabla;
	other.tabla = tabla;
	tabla = tabla_t;
	//
	votObjs.swap(other.votObjs);
}


void L_TransfSimil2DCeldaHash::voto(L_Calce &cal, L_List<L_Calce> *cal_pool) // Vota por los (i,j,k,z) compatibles con los (x0,y0,ang,sigma0) dados
{
	L_PUSH_EXECUTING_FN("L_TransfSimil2DCeldaHash::voto");
	int i0, i1, j0, j1, k0, k1, z;
	double ifrac, jfrac, kfrac, zfrac;
	double e,cosdt,sindt;

	votObjs.addVoteToFrequencies(cal.dRef->nobj);

	if (tabla==NULL)
	{
#ifdef L_BASIC_DEBUG
		if (nTabla_param<=0)
			L_hard_shutdown("L_TransfSimil2DCeldaHash::voto : non-positive size allocation");
#endif
		tabla=new L_TransfSimil2DCeldaLista[nTabla_param];
		nTabla=nTabla_param;
	}

	// zfrac no depende de nada
	// kfrac no depende de nada
	// ifrac depende de zfrac
	// jfrac depende de zfrac

	throw_L_ArgException_if(cal.dRef == NULL || cal.dPru == NULL, "L_TransfSimil2DCeldaHash::voto() : calce indefinido");

	zfrac = ( log(cal.dRef->sigma0)-log(cal.dPru->sigma0) )/L_LOG_2 / doctBin;  //  log(x) / log(2)
	kfrac = (cal.dRef->ang-cal.dPru->ang)/dangBin; // ver por cuales k votar
	cosdt = cos(cal.dRef->ang-cal.dPru->ang);
	sindt = sin(cal.dRef->ang-cal.dPru->ang);
#define L_ENTERO_V(x) ((int)((x)+0.5+10000)-10000)

	z = L_ENTERO_V(zfrac-0.5);
	e = pow(2.0,z);
	ifrac = ( e*cal.dPru->x0-cal.dRef->x0*cosdt+cal.dRef->y0*sindt )/(dxBin*cal.dRef->imorig_lx*e);
	jfrac = ( e*cal.dPru->y0-cal.dRef->y0*cosdt-cal.dRef->x0*sindt )/(dyBin*cal.dRef->imorig_ly*e);
	i0 = L_ENTERO_V(ifrac-0.5);
	i1 = L_ENTERO_V(ifrac+0.5);
	j0 = L_ENTERO_V(jfrac-0.5);
	j1 = L_ENTERO_V(jfrac+0.5);
	k0 = L_ENTERO_V(kfrac-0.5);
	k1 = L_ENTERO_V(kfrac+0.5);
	// 8 votos por los i,j,k,z enteros mas cercanos
	minivoto(i0, j0, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j0, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j1, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j1, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j0, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j0, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j1, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j1, k1, z, cal.dRef->nobj, cal, cal_pool);

	z = L_ENTERO_V(zfrac+0.5);
	e = pow(2.0,z);
	ifrac = ( e*cal.dPru->x0-cal.dRef->x0*cosdt+cal.dRef->y0*sindt )/(dxBin*cal.dRef->imorig_lx*e);
	jfrac = ( e*cal.dPru->y0-cal.dRef->y0*cosdt-cal.dRef->x0*sindt )/(dyBin*cal.dRef->imorig_ly*e);
	i0 = L_ENTERO_V(ifrac-0.5);
	i1 = L_ENTERO_V(ifrac+0.5);
	j0 = L_ENTERO_V(jfrac-0.5);
	j1 = L_ENTERO_V(jfrac+0.5);
	k0 = L_ENTERO_V(kfrac-0.5);
	k1 = L_ENTERO_V(kfrac+0.5);
	// 8 votos por los i,j,k,z enteros mas cercanos
	minivoto(i0, j0, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j0, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j1, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j1, k0, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j0, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j0, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i0, j1, k1, z, cal.dRef->nobj, cal, cal_pool);
	minivoto(i1, j1, k1, z, cal.dRef->nobj, cal, cal_pool);
#undef L_ENTERO_V
	L_POP_EXECUTING_FN("L_TransfSimil2DCeldaHash::voto");
}





void L_VerifPixelsTransfRGB::verifListaTransf(L_ImageRGBDouble &imPru, L_TransfPuntInt2DListaPtr &trL)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfRGB::verifListaTransf");
	L_TransfPuntInt2DNodoPtr *tr;
	for (tr=trL.root; tr!=NULL; tr=tr->sig)
	{
		if (!verifTransf(imPru, *tr->c))
		{
			delete tr->c;
			tr->c=NULL;
		}
	}
	trL.n=0;
	for (trL.pult=&trL.root; *trL.pult!=NULL;)
	{
		if ((*trL.pult)->c==NULL)
		{
			tr=*trL.pult;
			*trL.pult=(*trL.pult)->sig;
			delete(tr);
		}
		else
		{
			trL.pult=&(*trL.pult)->sig;
			trL.n++;
		}
	}
	L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::verifListaTransf");
}

void L_VerifPixelsTransfRGB::liberaPixeles()
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfRGB::liberaPixeles");
	L_ImagenRGBNodo *ptr;
	for (ptr=imL.root; ptr!=NULL; ptr=ptr->sig)
		ptr->c.destroy();
	L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::liberaPixeles");
}




void L_Debug_filtrosTransf::agregaFiltro(const char *name)
{
	L_PUSH_EXECUTING_FN("L_Debug_filtrosTransf::agregaFiltro");
	if (nFiltros >= L_DEFT_NFIL)
	{
		L_POP_EXECUTING_FN("L_Debug_filtrosTransf::agregaFiltro");
		return;
	}
	if (strlen(name)>L_DEFT_NCAR)
	{
		L_POP_EXECUTING_FN("L_Debug_filtrosTransf::agregaFiltro");
		return;
	}
	strcpy(filtros[nFiltros], name);
	nFiltros++;
	L_POP_EXECUTING_FN("L_Debug_filtrosTransf::agregaFiltro");
}

void L_Debug_filtrosTransf::marcaEliminacion_(L_TransfPuntInt2D *candTr, L_TransfSimil2DCelda *candCel, const char *name, double razon)
{
	L_PUSH_EXECUTING_FN("L_Debug_filtrosTransf::marcaEliminacion_");
	int i;
	for (i=0; i<nFiltros; i++)
	{
		if (strcmp(name,filtros[i])!=0)
			continue;
		switch(selPunt)
		{
		case nVotos:
			if (candTr!=NULL)
				insertaElimDebug(i, candTr->calL.size(), candTr->calL.size(), razon);
			else if (candCel!=NULL)
				insertaElimDebug(i, candCel->calL.size(), candCel->calL.size(), razon);
			break;
		case probTransf:
			throw L_NonImplementedException();
			break;
		case probCelda:
			throw L_NonImplementedException();
			break;
		}
		break;
	}
	L_POP_EXECUTING_FN("L_Debug_filtrosTransf::marcaEliminacion_");
}

void L_Debug_filtrosTransf::insertaElimDebug(int filtro, double puntaje, double caract, double razon)
{
	L_PUSH_EXECUTING_FN("L_Debug_filtrosTransf::insertaElimDebug");
	int i, j;
	nTransfTot++;
	if (nTransf >= L_DEFT_NPUN)
	{
		L_POP_EXECUTING_FN("L_Debug_filtrosTransf::insertaElimDebug");
		return;
	}
	for (i=0; i<nTransf; i++)
	{
		if (puntaje>mejoresPunt[i])
		{
			for(j=nTransf-1; j>=i; j--)
			{
				mejoresPunt[j+1]=mejoresPunt[j];
				mejoresCaract[j+1]=mejoresCaract[j];
				mejoresFiltros[j+1]=mejoresFiltros[j];
				mejoresRazon[j+1]=mejoresRazon[j];
			}
			mejoresPunt[i]=puntaje;
			mejoresCaract[i]=caract;
			mejoresFiltros[i]=filtro;
			mejoresRazon[i]=razon;
			nTransf++;
			L_POP_EXECUTING_FN("L_Debug_filtrosTransf::insertaElimDebug");
			return;
		}
	}
	mejoresPunt[nTransf]=puntaje;
	mejoresCaract[nTransf]=caract;
	mejoresFiltros[nTransf]=filtro;
	mejoresRazon[nTransf]=razon;
	nTransf++;
	L_POP_EXECUTING_FN("L_Debug_filtrosTransf::insertaElimDebug");
	return;
}

void L_Debug_filtrosTransf::fusiona(const L_Debug_filtrosTransf &d1, const L_Debug_filtrosTransf &d2)
{
	int i, u, v;
	bool agr1, agr2;

	nFiltros = d1.nFiltros + d2.nFiltros;
	for (i=0; i<d1.nFiltros; i++)
		strcpy(filtros[i],d1.filtros[i]);
	for (i=0; i<d2.nFiltros; i++)
		strcpy(filtros[d1.nFiltros+i],d2.filtros[i]);

	u=0;
	v=0;
	for (i=0; i<L_DEFT_NPUN; i++)
	{
		if (u < d1.nTransf)
		{
			if (v < d2.nTransf)
			{
				if (d1.mejoresPunt[u] > d2.mejoresPunt[v])
				{
					agr1=true;
					agr2=false;
				}
				else
				{
					agr1=false;
					agr2=true;
				}
			}
			else
			{
				agr1 = true;
				agr2 = false;
			}
		}
		else
		{
			if (v < d2.nTransf)
			{
				agr1 = false;
				agr2 = true;
			}
			else
			{
				agr1 = false;
				agr2 = false;
			}
		}
		if (agr1)
		{
				mejoresPunt[i] = d1.mejoresPunt[u];
				mejoresCaract[i] = d1.mejoresCaract[u];
				mejoresFiltros[i] = d1.mejoresFiltros[u];
				mejoresRazon[i] = d1.mejoresRazon[u];
				u++;
		}
		if (agr2)
		{
			mejoresPunt[i] = d2.mejoresPunt[v]+d1.nFiltros;
			mejoresCaract[i] = d2.mejoresCaract[v]+d1.nFiltros;
			mejoresFiltros[i] = d2.mejoresFiltros[v]+d1.nFiltros;
			mejoresRazon[i] = d2.mejoresRazon[v]+d1.nFiltros;
			v++;
		}
	}
	nTransf = u+v;
	nTransfTot = d1.nTransfTot + d2.nTransfTot;
	activo = d1.activo || d2.activo;
	selPunt = d1.selPunt; // Esto es arbitrario y poco importante
}

void L_Debug_filtrosTransf::imprimeReporte()
{
	L_PUSH_EXECUTING_FN("L_Debug_filtrosTransf::imprimeReporte");
	if (!activo)
	{
		L_POP_EXECUTING_FN("L_Debug_filtrosTransf::imprimeReporte");
		return;
	}
	int i;
	const char *caract;
	switch(selPunt)
	{
		case nVotos:
			caract="nVotos";
			break;
		case probTransf:
			caract="probTr";
			break;
		case probCelda:
			caract="probCel";
			break;
		default:
			caract="undef";
			break;
		}
	printf("Eliminaciones totales: %d\n", nTransfTot);
	for (i=0; i<nTransf; i++)
		printf(">%s: %.3f  filtro: %s: %.3f\n", caract, mejoresCaract[i], filtros[mejoresFiltros[i]], mejoresRazon[i]);
	L_POP_EXECUTING_FN("L_Debug_filtrosTransf::imprimeReporte");
}

bool L_VerifPixelsTransfBN::agregaImagenEntrenRobaObjetoDe(L_ImageGrayDouble &imRef)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenRobaObjetoDe");
	ultIm=imL.insertaRobaImagen(imRef);
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenRobaObjetoDe");
	return true;
}

bool L_VerifPixelsTransfBN::agregaImagenEntrenCopiaObjetoDe(L_ImageGrayDouble &imRef)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenCopiaObjetoDe");
	L_ImageGrayDouble other;
	other=imRef;
	ultIm=imL.insertaRobaImagen(other);
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenCopiaObjetoDe");
	return true;
}

bool L_VerifPixelsTransfBN::agregaImagenEntrenSinCopia(L_ImageGrayDouble &imRef)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenCopiaObjetoDe");
	ultIm=imL.insertaPtrImagen(imRef);
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::agregaImagenEntrenCopiaObjetoDe");
	return true;
}

double L_VerifPixelsTransfBN::calcCorrPix(const L_ImageGrayDouble &imPru, L_TransfPuntInt2D &tr)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::calcCorrPix");
	if (tr.calL.root==NULL)
	{
		L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::verifTransf");
		return false;
	}
	const L_ImageGrayDouble &imRef=*tr.calL.root->c.dRef->imorig;
	int i, j;
	int u, v;
	double lum11=0, lum12=0, lum22=0;
	double lum1=0, lum2=0;
	double a, b;
	double v11, v12, v22;
	double r;
	int nlum=0;
	L_PuntInt pRef, pPru;

	if (imRef.data()==NULL)
		throw L_ArgException();

	for (j=0; j<imRef.ly; j++)
	{
		for (i=0; i<imRef.lx; i++)
		{
			pRef.x0=i;
			pRef.y0=j;
			tr.proyeccionDe(pRef,pPru);
			u=(int)(pPru.x0+0.5);
			v=(int)(pPru.y0+0.5);
			if (u<0 || v<0 || u>=imPru.lx || v>=imPru.ly)
				continue;
			a=imRef.pix(i,j);
			b=imPru.pix(u,v);
			lum1+=a;
			lum2+=b;
			lum11+=a*a;
			lum12+=a*b;
			lum22+=b*b;
			nlum++;
		}
	}
	if (nlum==0)
		return 0; // Esto no deberia pasar...
	lum1/=nlum;
	lum2/=nlum;
	lum11/=nlum;
	lum12/=nlum;
	lum22/=nlum;
	v11=lum11 - lum1*lum1;
	v12=lum12 - lum1*lum2;
	v22=lum22 - lum2*lum2;

	r=v12/sqrt(v11*v22);
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::calcCorrPix");
	return r;
}

bool L_VerifPixelsTransfBN::verifTransf(const L_ImageGrayDouble &imPru, L_TransfPuntInt2D &tr, double *val)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::verifTransf");
	double corr = calcCorrPix(imPru, tr);
	if (val != NULL)
		*val = corr;
	if (corr >= umbCorr)
	{
		L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::verifTransf");
		return true;
	}
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::verifTransf");
	return false;
}

void L_VerifPixelsTransfBN::verifListaTransf(L_TransfPuntInt2DListaPtr &trL)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::verifListaTransf");
	L_TransfPuntInt2DNodoPtr *tr;
	const L_ImageGrayDouble *imPru=NULL;
	if (trL.size()==0)
		return;
	debFiltr.resetear();
	imPru=trL.root->c->calL.root->c.dPru->imorig;
	// Eliminar transaformaciones con poca correlacion de pixeles
	for (tr=trL.root; tr!=NULL; tr=tr->sig)
	{
		if (!verifTransf(*imPru, *tr->c, &tr->c->car.corrPix))
		{
			debFiltr.marcaEliminacion(tr->c, NULL, "corrPix", tr->c->car.corrPix);
			tr->c->calL.clear();
			delete tr->c;
			tr->c=NULL;
		}
	}
	// Rearmar la lista de transformaciones
	trL.n=0;
	for (trL.pult=&trL.root; *trL.pult!=NULL;)
	{
		if ((*trL.pult)->c==NULL)
		{
			tr=*trL.pult;
			*trL.pult=(*trL.pult)->sig;
			delete(tr);
		}
		else
		{
			trL.pult=&(*trL.pult)->sig;
			trL.n++;
		}
	}
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::verifListaTransf");
}

void L_VerifPixelsTransfBN::liberaPixeles()
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfBN::liberaPixeles");
	L_ImagenBNNodo *ptr;
	for (ptr=imL.root; ptr!=NULL; ptr=ptr->sig)
		ptr->c.destroy();
	L_POP_EXECUTING_FN("L_VerifPixelsTransfBN::liberaPixeles");
}


bool L_VerifPixelsTransfRGB::agregaImagenEntrenRobaObjetoDe(L_ImageRGBDouble &imRef)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfRGB::agregaImagenEntrenRobaObjetoDe");
	ultIm=&imL.push_back_swapping_ret(imRef)->c;
	L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::agregaImagenEntrenRobaObjetoDe");
	return true;
}

bool L_VerifPixelsTransfRGB::agregaImagenEntrenCopiaObjetoDe(L_ImageRGBDouble &imRef)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfRGB::agregaImagenEntrenCopiaObjetoDe");
	L_ImageRGBDouble other;
	other=imRef;
	ultIm=&imL.push_back_swapping_ret(other)->c;
	L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::agregaImagenEntrenCopiaObjetoDe");
	return true;
}

bool L_VerifPixelsTransfRGB::verifTransf(L_ImageRGBDouble &imPru, L_TransfPuntInt2D &tr)
{
	L_PUSH_EXECUTING_FN("L_VerifPixelsTransfRGB::verifTransf");
	if (tr.calL.root==NULL)
	{
		L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::verifTransf");
		return false;
	}
	const L_ImageRGBDouble &imRef=imPru;//*tr.calL.root->dRef->imorig; // REVISAR ESTA LINEA -> ERROR
	int i, j;
	int u, v;
	double lum11=0, lum12=0, lum22=0;
	double lum1R=0, lum2R=0;
	double lum1G=0, lum2G=0;
	double lum1B=0, lum2B=0;
	double a1, b1, a2, b2, a3, b3;
	double v11, v12, v22;
	double r;
	int nlum=0;
	L_PuntInt pRef, pPru;

	if (imRef.data()==NULL)
		throw L_ArgException();

	for (j=0; j<imRef.ly; j++)
	{
		for (i=0; i<imRef.lx; i++)
		{
			pRef.x0=i;
			pRef.y0=j;
			tr.proyeccionDe(pRef,pPru);
			u=(int)(pPru.x0+0.5);
			v=(int)(pPru.y0+0.5);
			if (u<0 || v<0 || u>=imPru.lx || v>=imPru.ly)
				continue;
			a1=imRef.pix(i,j,0);
			b1=imPru.pix(u,v,0);
			a2=imRef.pix(i,j,1);
			b2=imPru.pix(u,v,1);
			a3=imRef.pix(i,j,2);
			b3=imPru.pix(u,v,2);

			lum1R+=a1;
			lum1G+=a2;
			lum1B+=a3;
			lum2R+=b1;
			lum2G+=b2;
			lum2B+=b3;

			lum11+=a1*a1+a2*a2+a3*a3;
			lum12+=a1*b1+a2*b2*a3*b3;
			lum22+=b1*b1+b2*b2+b3*b3;
			nlum++;
		}
	}
	lum1R/=nlum;
	lum1G/=nlum;
	lum1B/=nlum;
	lum2R/=nlum;
	lum2G/=nlum;
	lum2B/=nlum;
	lum11/=nlum;
	lum12/=nlum;
	lum22/=nlum;
	v11=(lum11 - (lum1R*lum1R+lum1G*lum1G+lum1B*lum1B));
	v12=(lum12 - (lum1R*lum2R+lum1G*lum2G+lum1B*lum2B));
	v22=(lum22 - (lum2R*lum2R+lum2G*lum2G+lum2B*lum2B));

	r=v12/sqrt(v11*v22);
	if (r >= umbCorr)
	{
		L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::verifTransf");
		return true;
	}
	L_POP_EXECUTING_FN("L_VerifPixelsTransfRGB::verifTransf");
	return false;
}



bool L_GenTransfAfinHough::genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::genTransf");
	int i;
	L_CalceNodo *cal;
	L_TransfSimil2DCelda *bin;
	L_TransfSimil2DCeldaNodo *binNodo;
	L_TransfAfinPI2D tr;
	double lx, ly;
	double numDouble;
	double t1, t2, t3, t4;

	if (mostrarTiemposInterno)
	{
		t1 = L_TIME();
		printf("tiempos genTransf:  ");
	}
	debFiltr.resetear();

	ransacCoefRecalcular(); // Se debe llamar siempre por si se han modificado los parametros en forma externa

	for (cal=calL.root; cal!=NULL; cal=cal->sig)
		hash.voto(cal->c, cal_pool); // La lista calL no es alterada por este proceso

	if (mostrarTiemposInterno)
		t2 = L_TIME();

	// ACA : ordenar las celdas antes de procesarlas
	// Opcion de procesar una cantidad de celdas set

	L_Array<L_TransfSimil2DCeldaPtr> arrBin;

	if (hash.tabla!=NULL)
	{
		for (i=0; i<hash.nTabla; i++)
		{
			for (binNodo=hash.tabla[i].root; binNodo!=NULL; binNodo=binNodo->sig)
			{
				if ((int)binNodo->c.calL.size() >= nVotMinCelda)
				{
					arrBin.resize(arrBin.size()+1);
					arrBin[arrBin.size()-1].c = &binNodo->c;
				}
			}
		}
	}

	if (arrBin.size() > 0)
	{
		std::sort(&arrBin[0], &arrBin[0]+arrBin.size());

		if (nMaxCeldasRev != -1 && (int)arrBin.size() > nMaxCeldasRev)
			arrBin.resize(nMaxCeldasRev);

		// Revisar celdas
		for (int ic = 0; ic < (int)arrBin.size(); ic++)
		{
			bin = arrBin[ic].c;

			if (maxVot4D == true && ( // Se revisan sólo 4 direcciones "puras"
				(int)bin->calL.size() < hash.getVotos(bin->i+1, bin->j, bin->k, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i-1, bin->j, bin->k, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j+1, bin->k, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j-1, bin->k, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j, bin->k+1, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j, bin->k-1, bin->z, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j, bin->k, bin->z+1, bin->nobj) ||
				(int)bin->calL.size() < hash.getVotos(bin->i, bin->j, bin->k, bin->z-1, bin->nobj) )
				)
				goto RECHAZA_TRANSF;

			if (bin->calL.size() > (int)calL.size())
				bin->calL.eliminaRepetidos();

			tr.car.numVotosCelda = bin->calL.size();

			if (corrLinDatos.verif)
			{
				if (!verifCorrLinDatos(bin->calL, &numDouble))
				{
					debFiltr.marcaEliminacion(NULL, bin, "corrLinPuntos", numDouble);
					goto RECHAZA_TRANSF; // Este es un goto BIEN usado.
				}
				tr.car.corrLin = numDouble;
			}

			if (probCelda.verif)
			{
				numDouble=calcProbCelda(*bin, calL);
				if (numDouble<probCelda.pUmb)
				{
					debFiltr.marcaEliminacion(NULL, bin, "probCelda", numDouble);
					goto RECHAZA_TRANSF;
				}
				if (numDouble > 1)
					printf("calcProbCelda() > 1 con tr.calL.size() = %ld ...\n", long(tr.calL.size()));
				tr.car.probCelda = numDouble;
				tr.car.punt = numDouble;
			}

			bin->calL.moveListTo(tr.calL);

			if (!tr.calcTransfMinCuad()) // Se generan los parametros de la la transformacion
				goto RECHAZA_TRANSF;

			// Aqui deben ir todas las verificaciones
			if (distorRomb.verif)
			{
				if (!verifDistorRomb(tr, &numDouble))
				{
					debFiltr.marcaEliminacion(&tr, NULL, "distorRomb", numDouble);
					goto RECHAZA_TRANSF;
				}
				tr.car.distRomb = numDouble;
			}
			if (topDown==true)
				agregaTopDown(tr, calL, cal_pool); // POR ACA ESTA EL ERROR

			if (probTransf.verif==true)
			{
				numDouble = calcProbTransf(tr, calL);
				if (numDouble < probTransf.pUmb)
				{
					debFiltr.marcaEliminacion(&tr, NULL, "probTransf", numDouble);
					goto RECHAZA_TRANSF;
				}
				if (numDouble > 1)
					printf("calcProbTransf() > 1 con tr.calL.size() = %ld ...\n", long(tr.calL.size()));
				tr.car.probTransf = numDouble;
				tr.car.punt = numDouble;
			}
			lx=tr.calL.root->c.dRef->imorig_lx;
			ly=tr.calL.root->c.dRef->imorig_ly;
			if (ransac.activo && tr.calL.size()>=ransac.votMin && tr.calcTransfRansac(ransac.dx*lx, ransac.dy*ly, fnNumRansac((int)tr.calL.size()), ransac.nInt))
			{
				tr.car.ransac = true;
				tr.car.punt++; // Esto causaba acumulacion del puntaje en iteraciones sucesivas
			}
			else
				tr.car.ransac = false;

			tr.car.numVotosTransf = tr.calL.size();
			tr.car.numVotosTransfFus = tr.calL.size();

			// Finalmente se acepta la transformacion
			traL.push_back_swapping(tr);

	RECHAZA_TRANSF: // Este es un uso bacan de goto
			if (cal_pool != NULL)
				tr.calL.moveListTo(*cal_pool);
			else
				tr.calL.clear();
			tr.car.punt = 0;
		}
	}

	// clean tabla de hashing
	if (hash.tabla!=NULL)
		for (i=0; i<hash.nTabla; i++)
			hash.tabla[i].clear();

	if (mostrarTiemposInterno)
		t3 = L_TIME();

	// Fusionar transformacions similares
	traL.sort();
	if (fusion.activo)
		while(fusiona(traL, cal_pool));
	if (mostrarTiemposInterno)
		t4 = L_TIME();
	if (mostrarTiemposInterno)
		printf("%f  %f  %f  seg\n", t2-t1, t3-t2, t4-t3);

	// hash.destroy();  La tabla de hashing se destruyo 2 llaves mas arriba
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::genTransf");
	return true;
}

int L_GenTransfAfinHough::fnNumRansac(int nvotos)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::fnNumRansac");
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::fnNumRansac");
	return (int)(0.5 + nvotos * (ransac.porcVotInf + ransac.coef*pow(ransac.factorExp,nvotos)) );
}

void L_GenTransfAfinHough::agregaTopDown(L_TransfAfinPI2D &tr, L_CalceLista &calL, L_List<L_Calce> *cal_pool)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::agregaTopDown");
	L_CalceNodo *cal;
	L_TransfAfinPI2D itr;
	L_PuntInt res;
	if (tr.calL.root==NULL) // No debería pasar
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::agregaTopDown");
		return;
	}
	itr.calcParametrosInversosDe(tr);
	for (cal=calL.root; cal!=NULL; cal=cal->sig)
	{
		if (cal->c.dRef->nobj != tr.calL.root->c.dRef->nobj)
			continue;
		itr.proyeccionDe(*cal->c.dPru, res);
		if ( fabs((res.x0-cal->c.dRef->x0)/(double)cal->c.dRef->imorig_lx) >  hash.dxBin/2 )
			continue;
		if ( fabs((res.y0-cal->c.dRef->y0)/(double)cal->c.dRef->imorig_ly) >  hash.dyBin/2 )
			continue;
		tr.calL.push_back(cal->c, cal_pool);
	};
	tr.calL.sort(1, &L_CalceLista::cmpCalceXY);
	tr.calcTransfMinCuad();
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::agregaTopDown");
	return;
}

double L_GenTransfAfinHough::calcProbCelda(L_TransfSimil2DCelda &bin, L_CalceLista &calL, int nVotCelda) // Calculo de probabilidad nuevo
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::calcProbCelda");
	double rx, ry;
	double ca, sa;
	struct L_probComp{double d; double l; double r; double s;} p;
	double pVotoEnCelda=0;
	double p_celdaAzar;
	double p_celdaCorr;
	double ang;
	double imPru_lx,imPru_ly,imRef_lx,imRef_ly;
	imPru_lx=bin.calL.root->c.dPru->imorig_lx;
	imPru_ly=bin.calL.root->c.dPru->imorig_ly;
	imRef_lx=bin.calL.root->c.dRef->imorig_lx;
	imRef_ly=bin.calL.root->c.dRef->imorig_ly;
	ang=bin.k*hash.dangBin;
	ca=cos(ang);
	sa=sin(ang);

	// Explicacion:
	// Para esta test, el universo de los calces es la cantidad total de calces tentativos entre las 2 imagenes
	// La probabilidad de que un calce al azar vote por una celda (i,j,k,z) cuando vota por la mas cercana es:
	// P(i,j,k,z,obj) = P(nobj)*P(k)*P(z)*P(i,j|k,z)
	// P(nobj) = votos_hacia_esta_imagen_de_referencia / votos_totales
	// P(k) = dangBin / 2*Pi
	// P(z) = 0.6*pow(0.25, abs(z)) // Para width de 1 octava en z
	// p(i,j|k,z) = dxBin*dyBin / ( (rx+1) * (ry+1) )
	// pVotoAzar = P(nobj) * P(k) * P(z) * P(i,j|k,z)
	// pCeldaAzar depende de pVotoAzar (binomial)
	// p(deteccion) = pCalCorr / (pCalCorr + pCeldaAzar)
	//
	// La probabilidad de que un calce al azar vote por una celda (i,j,k,z) cuando vota por las 16 mas cercanas es:
	// pCeldaAzar = S(u=-1..0; v=-1..0; w=-1..0; x=-1..0) P(i+u,j+v,k+w,z+x,nobj)
	// pCeldaAzar = S(u=-1..0; v=-1..0; w=-1..0) ( P(i+u,j+v,k+w,z-1,nobj) + P(i+u,j+v,k+w,z,nobj) ) // No depende del k particular
	// pCeldaAzar = S(u=-1..0; v=-1..0) ( 2*P(i+u,j+v,k,z-1,nobj) + 2*P(i+u,j+v,k,z,nobj) ) // No depende del (i,j) particular
	// pCeldaAzar = 8*P(i,j,k,z-1,nobj) + 8*P(i,j,k,z,nobj)
	// pCeldaAzar = 8*P(i,j,k,z-1,nobj) + 8*P(i,j,k,z,nobj)
	// pCeldaAzar = 8*P(nobj)*P(k)*P(z-1)*P(i,j|k,z-1) + 8*P(nobj)*P(k)*P(z)*P(i,j|k,z)

	if (probCelda.calce16votos==true)
	{
		int dz;
		for (dz=-1; dz<1; dz++)
		{
			rx=sqrt(imPru_lx*imPru_lx*ca*ca+imPru_ly*imPru_ly*sa*sa)/(pow(2.0,bin.z+dz)*imRef_lx);
			ry=sqrt(imPru_ly*imPru_ly*ca*ca+imPru_lx*imPru_lx*sa*sa)/(pow(2.0,bin.z+dz)*imRef_ly);

			p.d=hash.votObjs.votosDe(bin.nobj)/hash.votObjs.nVotTot;
			p.l=hash.dxBin*hash.dyBin/((rx+1)*(ry+1)); //hash.dxBin=1/4
			p.r=hash.dangBin/(2*M_PI);
			p.s=3.0/5*pow(1.0/4,abs(bin.z+dz));
			pVotoEnCelda += 8 * p.d * p.l * p.r * p.s;
		}
	}
	else
	{
		rx=sqrt(imPru_lx*imPru_lx*ca*ca+imPru_ly*imPru_ly*sa*sa)/(pow(2.0,bin.z)*imRef_lx);
		ry=sqrt(imPru_ly*imPru_ly*ca*ca+imPru_lx*imPru_lx*sa*sa)/(pow(2.0,bin.z)*imRef_ly);

		p.d=hash.votObjs.votosDe(bin.nobj)/hash.votObjs.nVotTot;
		p.l=hash.dxBin*hash.dyBin/((rx+1)*(ry+1)); //hash.dxBin=1/4
		p.r=hash.dangBin/(2*M_PI);
		p.s=3.0/5*pow(1.0/4,abs(bin.z));
		pVotoEnCelda += p.d * p.l * p.r * p.s;
	}

	if (nVotCelda==-1) // Si no se ha entregado externamente la cantidad de votos que hay que considerar
		nVotCelda=(int)bin.calL.size();


	if (nVotCelda > (int)calL.size())
	{
		printf("Error en calcProbCelda(): Numero de votos en celda > numero de calces tentativos totales\n");
		p_celdaAzar=1;
	}
	else if (pVotoEnCelda > 1)
	{
		printf("Error en calcProbCelda(): Probabilidad de que un calce al azar vote por la celda es > 1\n");
		p_celdaAzar=1;
	}
	else
	{
		// Como se vota por los 2 mas cercanos en (i,j,k,z), el espacio real que cubre cada bin es el doble en cada dimension
		try
		{
			p_celdaAzar=L_binom(calL.size(), nVotCelda, pVotoEnCelda);
		}
		catch (L_ArgException)
		{
			printf("Error en rango de probabilidad en formula binomial\n");
			p_celdaAzar=1;
		}
	}
	double pCalCorr=bin.calL.root->c.dPru->refGen->calces_probCorrecto;
	p_celdaCorr=pCalCorr/(pCalCorr+p_celdaAzar);

	L_POP_EXECUTING_FN("L_GenTransfAfinHough::calcProbCelda");
 	return p_celdaCorr;
}

double L_GenTransfAfinHough::calcProbTransf(L_TransfAfinPI2D &tr, L_CalceLista &calL)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::calcProbTransf");
	struct L_probComp{double d; double l; double r; double s;} p;
	double pVotoEnCelda;
	double p_transfAzar;
	double p_transfCorr;
	L_CalceNodo *cal;
	L_TransfAfinPI2D itr;
	L_PuntInt res;
	int nTot=0;

	if (tr.calL.root==NULL) // No debería pasar
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::calcProbTransf");
		return 0;
	}

	// Explicacion:
	// Para esta test, se considera que el universo de los calces son los calces contenidos en el area proyectada
	// La probabilidad de que un calce al azar haya quedado asociado a la transformación es:
	// P(calce al azar en la transformacion) = P(esta imagen) * P(esta posicion) * P(esta orientacion) * P(esta escala)
	// P(esta imagen) = votos_hacia_esta_imagen_de_referencia / votos_totales
	// P(esta posicion) = dxBin * dyBin
	// P(esta orientacion) = dangBin / 2*Pi
	// P(esta escala) : segun Lowe es = 0.5
	// P(esta escala) : segun yo es = 0.6*pow(0.25, abs(int(z-0.5))) + 0.6*pow(0.25, abs(int(z+0.5))), con z la diferencia de escala
	// p(transformacion debida al azar) depende de P(calce al azar en la transformacion) (binomial)
	// p(deteccion) = pCalCorr / (pCalCorr + P(transformacion debida al azar))
	// Las suposiciones anteriores consideran que todos los calces se originaron a partir de Hough, pero eso no es exacto
	//  *Top-down matching: addFrom calces que no deberian haber quedado en el bin
	//  *Se vota por las 16 celdas mas cercanas en Hough, no por la mas cercana

	itr.calcParametrosInversosDe(tr);
	// tr:  imagen de referencia -> imagen de test
	// itr: imagen de test -> imagen de referencia
	for (cal=calL.root; cal!=NULL; cal=cal->sig)
	{
		itr.proyeccionDe(*cal->c.dPru, res);
		if ( res.x0>=0 && res.x0<=cal->c.dRef->imorig_lx && res.y0>=0 && res.y0<=cal->c.dRef->imorig_ly)
			nTot++;
	}
	p.d=hash.votObjs.votosDe(tr.calL.root->c.dRef->nobj)/hash.votObjs.nVotTot;
	p.l=hash.dxBin*hash.dyBin;
	p.r=hash.dangBin/(2*M_PI);

	if (probTransf.corrPz==true)
	{
		double det, z_dec;
		int z1, z2;
		det=tr.m11*tr.m22-tr.m12*tr.m21; // Es la proporcion entre el area original y la proyectada (si la transf no se deforma)
		z_dec=log(det)/log(2.0) /2; // /log2: para pasarlo a log2(x). /2: para sacar el termino cuadratico (se ahorra una sqrt)
		z1 = (int) (z_dec - 0.5);
		z2 = (int) (z_dec + 0.5);
		p.s= 3.0/5*pow(1.0/4,abs(z1)) + 3.0/5*pow(1.0/4,abs(z2));
	}
	else
		p.s=0.5;

	pVotoEnCelda = p.d * p.l * p.r * p.s;

	if ((int)tr.calL.size() > nTot) // Hay un calce en la transformacion que se sale de la imagen de referencia al retroproyectarse, nada tan extrano
	{
		//printf("Error en calcProbTransf(): Numero de votos en la transformacion > numero de calces tentativos entre las areas proyectadas\n");
		nTot = (int)tr.calL.size();
	}
	if (pVotoEnCelda > 1) // Esto si es extrano
	{
		printf("Error en calcProbTransf(): Probabilidad de que un calce al azar vote por la transformacion es > 1\n");
		pVotoEnCelda=1;
	}

	try
	{
		p_transfAzar=L_binom(nTot, tr.calL.size(), pVotoEnCelda);
	}
	catch (L_ArgException)
	{
		printf("calcProbTransf(): Error en rango de probabilidad en calculo binomial\n");
		p_transfAzar=1;
	}
	double pCalCorr=calL.root->c.dPru->refGen->calces_probCorrecto;
	p_transfCorr=pCalCorr/(pCalCorr+p_transfAzar);
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::calcProbTransf");
	return p_transfCorr;
}

bool L_GenTransfAfinHough::verifDistorRomb(L_TransfAfinPI2D &tr, double *val)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::verifDistorRomb");
	double ladoA, ladoB, relLados, relArea;

	// La transformacion afin transformUsing:
	// (1,0) en (m11,m12) + (tx,ty)
	// (0,1) en (m21,m22) + (tx,ty)
	// (0,0) en (0,0)     + (tx,ty)

	ladoA=sqrt(tr.m11*tr.m11 + tr.m21*tr.m21);
	ladoB=sqrt(tr.m12*tr.m12 + tr.m22*tr.m22);
	if (ladoA>ladoB)
		relLados=ladoA/ladoB;
	else
		relLados=ladoB/ladoA;
	if (val!=NULL)
		*val=relLados;
	if (relLados>distorRomb.lados)
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifDistorRomb");
		return false;
	}
	relArea=(tr.m11*tr.m22-tr.m21*tr.m12)/(ladoA*ladoB);
	if (val!=NULL)
		*val=relArea;
	if (relArea<distorRomb.senAng)
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifDistorRomb");
		return false;
	}
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifDistorRomb");
	return true;
}

bool L_GenTransfAfinHough::fusiona(L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool)
{
	L_PUSH_EXECUTING_FN("bool L_GenTransfAfinHough::fusiona");
	L_TransfPuntInt2DNodoPtr *trP1, *trP2;
	L_TransfAfinPI2D *ta1, *ta2;
	L_TransfAfinPI2D *tr;
	double lx, ly;
	bool fusionados=false;
	int cmp=0;

	if (traL.root==NULL || traL.root->sig==NULL) // 1 transformacion o menos
	{
		L_POP_EXECUTING_FN("bool L_GenTransfAfinHough::fusiona");
		return false;
	}

	//for (trP1 = traL.root; trP1->sig!=NULL; trP1=trP1->sig)
	//	if (trP1->c->calL.verifNumElem() == false)
	//		system("echo error antes de fusionar >> err.txt");

	for (trP1=traL.root; trP1->sig!=NULL; trP1=trP1->sig)
	{
		for (trP2=trP1->sig; trP2!=NULL; trP2=trP2->sig)
		{
			ta1 = dynamic_cast<L_TransfAfinPI2D*>(trP1->c);
			ta2 = dynamic_cast<L_TransfAfinPI2D*>(trP2->c);

			if (ta1!=NULL && ta2!=NULL && verifFusion(*ta1, *ta2))
			{
				throw_L_ArgException_if(ta1->calL.root == NULL || ta2->calL.root == NULL, "L_GenTransfAfinHough::fusiona() : transformacion vacia no nula"); // Mejor que se caiga de una forma ordenada si pasa algo malo
				lx=trP1->c->calL.root->c.dRef->imorig_lx;
				ly=trP1->c->calL.root->c.dRef->imorig_ly;
				if (ta1->car.ransac && fusion.separarRansac)
				{
					if (ta2->car.ransac && fusion.separarRansac)
					{
						// ta1 y ta2 son ransac
						// calcular fusion y ver si es mejor que los padres, si no, clean el peor padre
						tr=new L_TransfAfinPI2D();
						//if (
						//	ta1->calL.copyListOn(tr->calL) == false ||
						//	ta2->calL.copyListOn(tr->calL) == false)
						//		system("echo error en medio de la fusion v1 >> err.txt");
						ta1->calL.copyListOn(tr->calL, cal_pool);
						ta2->calL.copyListOn(tr->calL, cal_pool);

						// En teoria deberian haber por lo menos 3 calces independientes
						tr->calL.sort(1, &L_CalceLista::cmpCalceXY, cal_pool); // Error: decia 2
						if (tr->calL.size()>ta1->calL.size() && tr->calL.size()>ta2->calL.size() &&
							tr->calcTransfMinCuad() && tr->calL.size() >= ransac.votMin &&
							tr->calcTransfRansac(ransac.dx*lx, ransac.dy*ly, fnNumRansac((int)tr->calL.size()), ransac.nInt) &&
							tr->calL.size()>ta1->calL.size() && tr->calL.size()>ta2->calL.size())
						{
							if (ta1->car.punt > ta2->car.punt || (ta1->car.punt==ta2->car.punt && ta1->calL.size() > ta2->calL.size()))
								tr->car.punt=ta1->car.punt;
							else
								tr->car.punt=ta2->car.punt;
							tr->car.ransac=true;
							tr->car.numVotosTransfFus = tr->calL.size();
							if (cal_pool != NULL)
							{
								ta1->calL.moveListTo(*cal_pool);
								ta2->calL.moveListTo(*cal_pool);
							}
							delete ta1;
							delete ta2;
							trP1->c=tr;
							trP2->c=NULL;
							fusionados=true;
						}
						else
						{
							delete tr;
							if (ta1->car.punt > ta2->car.punt || (ta1->car.punt==ta2->car.punt && ta1->calL.size() > ta2->calL.size()))
								cmp=1;
							else if (ta1->car.punt < ta2->car.punt)
								cmp=-1;
							else if (ta1->calL.size() > ta2->calL.size())
								cmp=1;
							else
								cmp=-1;
							if (cmp==1)
							{
								debFiltr.marcaEliminacion(ta2, NULL, "fusion_RANSAC", -2);
								if (cal_pool != NULL)
									ta2->calL.moveListTo(*cal_pool);
								delete ta2;
								trP2->c=NULL;
							}
							else
							{
								debFiltr.marcaEliminacion(ta1, NULL, "fusion_RANSAC", -2);
								if (cal_pool != NULL)
									ta1->calL.moveListTo(*cal_pool);
								delete ta1;
								trP1->c=NULL;
							}
						}
					}
					else
					{
						// ta1 es ransac solamente
						// clean la transformacion que no es ransac
						debFiltr.marcaEliminacion(ta2, NULL, "fusion_RANSAC", -2);
						if (cal_pool != NULL)
							ta2->calL.moveListTo(*cal_pool);
						delete ta2;
						trP2->c=NULL;
					}
				}//if (ta1->ransac)
				else
				{
					if (ta2->car.ransac && fusion.separarRansac)
					{
						// ta2 es ransac solamente
						// clean la transformacion que no es ransac
						debFiltr.marcaEliminacion(ta1, NULL, "fusion_RANSAC", -2);
						if (cal_pool != NULL)
							ta1->calL.moveListTo(*cal_pool);
						delete ta1;
						trP1->c=NULL;
					}
					else
					{
						// ninguno es ransac
						// calcular fusion y ver si es mejor que los padres, si no, clean el peor padre
						tr=new L_TransfAfinPI2D();
						//if (
						//	ta1->calL.copyListOn(tr->calL) == false ||
						//	ta2->calL.copyListOn(tr->calL) == false)
						//		system("echo error en medio de la fusion v2 >> err.txt");
						ta1->calL.copyListOn(tr->calL, cal_pool);
						ta2->calL.copyListOn(tr->calL, cal_pool);
						tr->calL.sort(1, &L_CalceLista::cmpCalceXY); // Error: decia 2
						if (tr->calL.size()>ta1->calL.size() && tr->calL.size()>ta2->calL.size() && tr->calcTransfMinCuad())
						{
							if (ta1->car.punt>ta2->car.punt)
								tr->car.punt=ta1->car.punt;
							else
								tr->car.punt=ta2->car.punt;
							if (cal_pool != NULL)
							{
								ta1->calL.moveListTo(*cal_pool);
								ta2->calL.moveListTo(*cal_pool);
							}
							delete ta1;
							delete ta2;
							trP1->c=tr;
							trP2->c=NULL;
							fusionados=true;
						}
						else
						{
							delete tr;
							if (ta1->car.punt > ta2->car.punt)
								cmp=1;
							else if (ta1->car.punt < ta2->car.punt)
								cmp=-1;
							else if (ta1->calL.size() > ta2->calL.size())
								cmp=1;
							else
								cmp=-1;
							if (cmp==1)
							{
								debFiltr.marcaEliminacion(ta2, NULL, "fusion_RANSAC", -1);
								if (cal_pool != NULL)
									ta2->calL.moveListTo(*cal_pool);
								delete ta2;
								trP2->c=NULL;
							}
							else
							{
								debFiltr.marcaEliminacion(ta1, NULL, "fusion_RANSAC", -1);
								if (cal_pool != NULL)
									ta1->calL.moveListTo(*cal_pool);
								delete ta1;
								trP1->c=NULL;
							}
						}
					}
				}//else (ta1->ransac)
			}//if (ta1!=NULL && ta2!=NULL && verifFusion(*ta1, *ta2))
		}// for (trP2)
	}//for (trP1)
	traL.n=0;
	for (traL.pult=&traL.root; *traL.pult!=NULL;)
	{
		if ((*traL.pult)->c==NULL)
		{
			// eliminar nodo
			trP1=*traL.pult;
			*traL.pult=(*traL.pult)->sig;
			delete(trP1);
		}
		else
		{
			// nada
			traL.pult=&(*traL.pult)->sig;
			traL.n++;
		}
	}
	L_POP_EXECUTING_FN("bool L_GenTransfAfinHough::fusiona");
	return fusionados;
}

bool L_GenTransfAfinHough::verifFusion(L_TransfAfinPI2D &tr1, L_TransfAfinPI2D &tr2)
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::verifFusion");
	double x1[4], y1[4], x2[4], y2[4];
	double lx, ly;
	double perim=0, dist=0, superp;
	int k, k2;

	if (tr1.calL.root->c.dPru->nobj != tr2.calL.root->c.dPru->nobj)
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifFusion");
		return false; // Corresponden a distintos objetos
	}

	lx=tr1.calL.root->c.dPru->imorig_lx;
	ly=tr1.calL.root->c.dPru->imorig_ly;

	x1[0]=tr1.tx; //TRx(tr1,0,0)
	y1[0]=tr1.ty; //TRy(tr1,0,0)
	x1[1]=tr1.m11*lx+tr1.tx; //TRx(tr1,lx,0);
	y1[1]=tr1.m21*lx+tr1.ty; //TRy(tr1,lx,0);
	x1[2]=tr1.m11*lx+tr1.m12*ly+tr1.tx; //TRx(tr1,lx,ly);
	y1[2]=tr1.m21*lx+tr1.m22*ly+tr1.ty; //TRy(tr1,lx,ly);
	x1[3]=tr1.m12*ly+tr1.tx; //TRx(tr1,0,ly);
	y1[3]=tr1.m22*ly+tr1.ty; //TRy(tr1,0,ly);

	x2[0]=tr2.tx; //TRx(tr2,0,0)
	y2[0]=tr2.ty; //TRy(tr2,0,0)
	x2[1]=tr2.m11*lx+tr2.tx; //TRx(tr2,lx,0);
	y2[1]=tr2.m21*lx+tr2.ty; //TRy(tr2,lx,0);
	x2[2]=tr2.m11*lx+tr2.m12*ly+tr2.tx; //TRx(tr2,lx,ly);
	y2[2]=tr2.m21*lx+tr2.m22*ly+tr2.ty; //TRy(tr2,lx,ly);
	x2[3]=tr2.m12*ly+tr2.tx; //TRx(tr2,0,ly);
	y2[3]=tr2.m22*ly+tr2.ty; //TRy(tr2,0,ly);

	for (k=0; k<4; k++)
	{
		if (k!=3) k2=k+1; else k2=0;
		perim+=sqrt((x1[k]-x1[k2])*(x1[k]-x1[k2])+(y1[k]-y1[k2])*(y1[k]-y1[k2]));
		perim+=sqrt((x2[k]-x2[k2])*(x2[k]-x2[k2])+(y2[k]-y2[k2])*(y2[k]-y2[k2]));
		dist+=sqrt((x1[k]-x2[k])*(x1[k]-x2[k])+(y1[k]-y2[k])*(y1[k]-y2[k]));
	}

	superp=1-(dist/perim);
	if (superp<fusion.superp)
	{
		L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifFusion");
		return false;
	}
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::verifFusion");
	return true;
}

void L_GenTransfAfinHough::swap(L_GenTransfAfinHough &other)
{
	other.L_GenTransfCalces::swap(*this);
	other.L_GenTransfAfinHough_POD::swap(*this);
	paramGTAH.swap(other.paramGTAH);
	hash.swap(other.hash);
}

L_ParamManagerLocal *L_GenTransfAfinHough::pideParams()
{
	L_PUSH_EXECUTING_FN("L_GenTransfAfinHough::pideParams");
	L_POP_EXECUTING_FN("L_GenTransfAfinHough::pideParams");
	return &paramGTAH;
}

L_GenTransfProyectivaRansac_1tr::L_GenTransfProyectivaRansac_1tr() : L_GenTransfCalces(genTrProyRansac_1tr), paramGTPR_1("GenTransfProyRansac", 3)
{
	radioError = 6;
	usar_4calces=true;
	considerarCam=false;
	minCalces=12;
	paramGTPR_1.addFrom("radioError", &radioError);
	paramGTPR_1.addFrom("usar_4calces", &usar_4calces);
	paramGTPR_1.addFrom("minCalces", &minCalces);
}

bool L_GenTransfProyectivaRansac_1tr::genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool)
{
	L_TransfProyectivaPI2D tr;
	int i;
	bool ret=false;
	tr.ransac_2cal=!usar_4calces;
	tr.camInter=cam;
	tr.distorsionarCoordImagenReferencia = considerarCam;
	for (i=0; i<3; i++)
	{
		if (cal_pool != NULL)
			tr.calL.moveListTo(*cal_pool);
		else
			tr.calL.clear();
		calL.copyListOn(tr.calL, cal_pool);
		ret = tr.calcTransfRansac(radioError, radioError, L_MAX(minCalces, (int)calL.size()/2), 30);
		if (ret == false)
			break;
		if (corrLinDatos.verif && verifCorrLinDatos(tr.calL) == false)
			continue; // Probar de nuevo, podria funcionar
		traL.push_back_swapping(tr);
	}
	return ret;
}

L_GenTransfCamParalelas_1tr::L_GenTransfCamParalelas_1tr() : L_GenTransfCalces(genTrCamParal_1tr), paramGTCP_1("GTCP", 3)
{
	errorVertAnt = 15;
	errorVertPost = 4;
	paramGTCP_1.addFrom("errorVertAnt", &errorVertAnt);
	paramGTCP_1.addFrom("errorVertPost", &errorVertPost);
}

bool L_GenTransfCamParalelas_1tr::genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool)
{
	L_CalceNodo **c, *ctmp;
	L_TransfAfinPI2D tr;
	if (calL.size() < 10)
		return false;
	calL.copyListOn(tr.calL);
	int celdas[12]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int i, iMax, valMax;
	int v;
	double dy;
	double dang;

	// Inicio eliminacion de nodo con condicion
	for (c=&tr.calL.root; *c!=NULL;)
	{
		dy = (*c)->c.dPru->y0 - (*c)->c.dRef->y0;
		dang = (*c)->c.dPru->ang - (*c)->c.dRef->ang;
		while (dang > M_PI)
			dang-=M_PI;
		while (dang < -M_PI)
			dang+=M_PI;
		if (
			L_ABS(dy) <= errorVertAnt
			&&
			(*c)->c.dPru->sigma0 / (*c)->c.dRef->sigma0 < 1.5
			&&
			(*c)->c.dPru->sigma0 / (*c)->c.dRef->sigma0 > 1.0/1.5
			&&
			L_ABS(dang) <= 30*M_PI/180
			)
		{
			v = 5 + (int)( dy * 5.0 / errorVertAnt );
			celdas[v]++;
			celdas[v+1]++;
			c=&(*c)->sig; // Avanzar en la lista
		}
		else
		{
			ctmp=*c;
			*c=(*c)->sig; // Desencadenar nodo
			if (cal_pool != NULL)
				cal_pool->push_back_node(ctmp);
			else
				delete ctmp; // clean nodo
			tr.calL.n--;
		}
	}
	tr.calL.pult=c; // Puntero al final de la lista
	// Fin eliminacion de nodo con condicion

	valMax = 0;
	for (i=0; i<10; i++)
	{
		if (celdas[i] > valMax)
		{
			valMax = celdas[i];
			iMax = i;
		}
	}

	// Inicio eliminacion de nodo con condicion
	for (c=&tr.calL.root; *c!=NULL;)
	{
		dy = (*c)->c.dPru->y0 - (*c)->c.dRef->y0;
		v = 5 + (int)( dy * 5.0 / errorVertAnt );
		if ( v == iMax || v+1 == iMax )
			c=&(*c)->sig; // Avanzar en la lista
		else
		{
			ctmp=*c;
			*c=(*c)->sig; // Desencadenar nodo
			if (cal_pool != NULL)
				cal_pool->push_back_node(ctmp);
			else
				delete ctmp; // clean nodo
			tr.calL.n--;
		}
	}
	tr.calL.pult=c; // Puntero al final de la lista
	// Fin eliminacion de nodo con condicion

	if ((int)tr.calL.size() > 4)
	{
		if (tr.calcTransfMinCuad() == true)
		{
			traL.push_back_swapping(tr);
			return true;
		}
	}


	return false;
}

bool L_RecPuntInt::agregaImagenEntren_interno(L_ImageGrayDouble &imRef, int modoAlmacenamiento, double s0, L_DescriptorLista *precalc)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::agregaImagenEntren_interno");
	ti.t1=L_TIME();
	gpi->nobj=nobj;
	switch(modoAlmacenamiento)
	{
	case 0: // Robar la imagen entregada
		corrPixBN.agregaImagenEntrenRobaObjetoDe(imRef);
		break;
	case 1: // Copiar la imagen entregada
		corrPixBN.agregaImagenEntrenCopiaObjetoDe(imRef);
		break;
	case 2: // Guardar referencia a la imagen entregada
		corrPixBN.agregaImagenEntrenSinCopia(imRef); // En este caso hay que tener cuidado de que la imagen pasada sobreviva lo suficiente
		break;
	}
	imsref.resize(imsref.size()+1);
	imsref[nobj++]=corrPixBN.ultIm;
	if (precalc==NULL)
	{
		gpi->calcDescrPuntInt(*corrPixBN.ultIm, s0);
		precalc=&gpi->desL;
	}
	else
	{
		L_DescriptorNodo *ptr;
		for (ptr=precalc->root; ptr!=NULL; ptr=ptr->sig)
		{
			ptr->c.imorig=corrPixBN.ultIm;
			ptr->c.nobj=gpi->nobj;
		}
	}
	precalc->moveListTo(desLRef);
	ti.t2=L_TIME();
	ti.gpiTime.push(ti.t2-ti.t1);
	L_POP_EXECUTING_FN("L_RecPuntInt::agregaImagenEntren_interno");
	return true;
}

bool L_RecPuntInt::agregaImagenesEntren_usa_desL(L_DescriptorLista &desL, L_ImageGrayDouble *imRef, bool copiarImagen)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::agregaImagenesEntren_usa_desL");
	ti.t1=L_TIME();
	conoceImagenesEntrenamiento = false;
	if (corrPixBN.verif == true)
	{
		if (imRef == NULL)
		{
			printf("Error: No se puede activar corrPixBN.verif si se entrena con agregaImagenesEntren_usa_desL(desL)\n");
			printf("corrPixBN.verif desactivado\n");
			corrPixBN.verif=false;
		}
	}

	if (copiarImagen)
	{
		corrPixBN.agregaImagenEntrenCopiaObjetoDe(*imRef);
		imsref.resize(imsref.size()+1);
		desL.moveListTo(desLRef);
		desLRef.fijaNumObj(nobj);
		imsref[nobj++]=corrPixBN.ultIm;
	}
	else // Aca la imagen que se le pasa debe preserarse sin cambios
	{
		imsref.resize(imsref.size()+1);
		desL.moveListTo(desLRef);
		desLRef.fijaNumObj(nobj);
		imsref[nobj++]=imRef;
	}

	ti.t2=L_TIME();
	ti.gpiTime.push(ti.t2-ti.t1);
	L_POP_EXECUTING_FN("L_RecPuntInt::agregaImagenesEntren_usa_desL");
	return true;
}

bool L_RecPuntInt::finalizaEntren()
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::finalizaEntren");
	ti.t1=L_TIME();
	if (filtD!=NULL)
		filtD->filtrarNodos(desLRef);
	gbd->createFrom(desLRef);
	ti.t2=L_TIME();
	ti.gbdTime.push(ti.t2-ti.t1);
	L_POP_EXECUTING_FN("L_RecPuntInt::finalizaEntren");
	return true;
}

bool L_RecPuntInt::procesaImagenesPrueba_arch(L_ImageGrayDouble &imPru, double s0, const char *dibCalcesDebug)
{
	L_ImageRGBUchar imCal;
	bool ret;
	ret = procesaImagenesPrueba(imPru, s0, &imCal);
	imCal.saveImage(dibCalcesDebug);
	return ret;
}

bool L_RecPuntInt::procesaImagenesPrueba(L_ImageGrayDouble &imPru, double s0, L_ImageRGBUchar *dibCalcesDebug)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::procesaImagenesPrueba(imPru,s0)");
	bool ret;

	ti.t1=L_TIME();
	calL.clear();
	gpi->nobj=-1;
	gpi->desL.clear();
	gpi->calcDescrPuntInt(imPru, s0);
	ti.t2=L_TIME();
	ti.gpiTime.push(ti.t2-ti.t1);

	desLPru.clear();
	gpi->desL.moveListTo(desLPru);

	ret = procesaImagenesPrueba_interno(&imPru, dibCalcesDebug); // Mide ti.gtcTime

	L_POP_EXECUTING_FN("L_RecPuntInt::procesaImagenesPrueba(imPru,s0)");
	return ret;
}

bool L_RecPuntInt::procesaImagenesPrueba_usa_desL_arch(const L_ImageGrayDouble *imPru, const char *dibCalcesDebug)
{
	bool ret;
	if (dibCalcesDebug != NULL)
	{
		L_ImageRGBUchar imCal;
		ret = procesaImagenesPrueba_interno(imPru, &imCal);
		imCal.saveImage(dibCalcesDebug);
	}
	else
		ret = procesaImagenesPrueba_interno(imPru);
	return ret;
}

bool L_RecPuntInt::procesaImagenesPrueba_usa_desL(L_DescriptorLista &desL, const L_ImageGrayDouble *imPru, const char *dibCalcesDebug)
{
	this->desLPru.clear();
	desL.moveListTo(this->desLPru);
	bool r=procesaImagenesPrueba_usa_desL_arch(imPru, dibCalcesDebug);
	return r;
} // Si se llama con imPru=NULL no se puede usar dibujaTransformaciones()


// La funcion que hace el trabajo
bool L_RecPuntInt::procesaImagenesPrueba_interno(const L_ImageGrayDouble *imPru, L_ImageRGBUchar *dibCalcesDebug)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::procesaImagenesPrueba()");
	bool ret;

	impru = imPru;

	if (cuentaDistrSigmas)
	{
		distrS.AgregaLista(desLPru);
		distrS.print(cuentaDistrS_imprDiv);
	}

	ti.t1=L_TIME();
	if (filtD!=NULL)
		filtD->filtrarNodos(desLPru);

	if (gpi->calces_numCalcesPorDescr==1)
		calL.genCalces(*gbd, desLPru, deltaAngCalces.apply, deltaAngCalces.deltaMax, &cal_pool);
	else
		calL.genCalces_nMasCercanos(*gbd, desLPru, deltaAngCalces.apply, deltaAngCalces.deltaMax, &cal_pool);

	//calL.revisaCoherencia();

	ti.t2=L_TIME();
	caractRec.nCal = (int)calL.size();
	ti.calTime.push(ti.t2-ti.t1);
	if (dibCalcesDebug != NULL)
		calL.dibuja(*dibCalcesDebug);

	ti.t1=L_TIME();
	traL.destroyList(&cal_pool);

	ret=gtc->genTransf(calL, traL, &cal_pool);

	if (RSL.apply)
		traL.aplicaRSL(RSL.corrige, RSL.nVecinos, RSL.minCorresp, RSL.verifEscalaAngs, RSL.varEscala, RSL.varAng);
	if (corrPixBN.verif)
		corrPixBN.verifListaTransf(traL); // Esto puede eliminar transformaciones
	if (cuenta_calces_1_a_n)
		traL.cuenta_calces_1_a_n();
	if (cuentaCerros.apply)
		traL.aplicaCuentaCerros(cuentaCerros.corrige, cuentaCerros.minDistCerros, cuentaCerros.maxDistCerros, cuentaCerros.minCerrosUtil, cuentaCerros.parecidoCerrosMin, cuentaCerros.nLineasMax, cuentaCerros.porcCerrosCorrectosMin);

	traL.sort();

	if (gtc->debFiltr.activo || corrPixBN.debFiltr.activo)
	{
		L_Debug_filtrosTransf deb;
		deb.fusiona(gtc->debFiltr, corrPixBN.debFiltr);
		deb.imprimeReporte();
	}

	ti.t2=L_TIME();
	caractRec.nTr = (int)traL.size();
	ti.gtcTime.push(ti.t2-ti.t1);

	L_POP_EXECUTING_FN("L_RecPuntInt::procesaImagenesPrueba()");
	return ret;
}

bool L_RecPuntInt::procesaReferenciaPruebaExt(L_ImageGrayDouble &imRef, L_ImageGrayDouble &imPru, double s0, L_ImageRGBUchar *dibCalcesDebug)
{
	if (agregaImagenEntren_interno(imRef,2,s0, NULL) == false)
		return false;
	finalizaEntren();
	return procesaImagenesPrueba(imPru,s0,dibCalcesDebug);
}

void L_RecPuntInt::reseteaEntren()
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::reseteaEntren");
	reseteaPrueba();
	this->nobj=0;
	gpi->desL.clear();
	desLRef.clear();
	desLPru.clear();
	if (gbd!=NULL)
		gbd->destroy();
	traL.destroyList(&cal_pool);
	impru=NULL;
	imsref.resize(0);
	corrPixBN.liberaPixeles();
	conoceImagenesEntrenamiento = true;
	L_POP_EXECUTING_FN("L_RecPuntInt::reseteaEntren");
}

bool L_RecPuntInt::dibujaTransformaciones(L_ImageRGBUchar &imSal, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica, int transfEspecifica)
{
	L_ShapeArray lins;
	L_Array <const L_ImageGrayDouble *> imarr;

	if (conoceImagenesEntrenamiento == false || impru == NULL) // El reconocedor no tiene acceso a todas las imagenes de entrenamiento y test
	{
		imSal.reallocate(1,1);
		printf("Advertencia: L_RecPuntInt::dibujaTransformaciones() no conoce las imagenes necesarias\n");
		return false;
	}

	traL.dibujaTransformaciones(lins, calces1_o_caja0,  maxImRef, maxTransfPorImagen, maxTransfTotales, imRefEspecifica, transfEspecifica);
	generaArregloImagenesParaCalces(imarr);
	imSal.genDrawingMatches(imarr, lins);
	return true;
}

void L_RecPuntInt::generaArregloImagenesParaCalces(L_Array<const L_ImageGrayDouble *> &imarr) const
{
	int i;
	imarr.resize(imsref.size()+1);;
	imarr[0]=impru;
	for (i=0; i<(int)imsref.size(); i++)
		imarr[i+1]=imsref[i];
}

bool L_RecPuntInt::dibujaCalcesIniciales(const char *nomarch)
{
	L_ImageRGBUchar im;
	dibujaCalcesIniciales(im);
#if __COMPAT_ATLIMAGE__
	return im.writeImageATL(nomarch)!=0;
#else
	return im.writeBMP(nomarch)!=0;
#endif
}

bool L_RecPuntInt::dibujaCalcesIniciales(L_ImageRGBUchar &imSal)
{
	L_ShapeArray lins;
	L_Array <const L_ImageGrayDouble *> imarr;

	if (conoceImagenesEntrenamiento == false || impru == NULL) // El reconocedor no tiene acceso a todas las imagenes de entrenamiento y test
	{
		imSal.reallocate(1,1);
		printf("Advertencia: L_RecPuntInt::dibujaTransformaciones() no conoce las imagenes necesarias\n");
		return false;
	}
	calL.genLineas(lins);
	generaArregloImagenesParaCalces(imarr);
	imSal.genDrawingMatches(imarr, lins);
	return true;
}

bool L_RecPuntInt::dibujaTransformaciones(const char *nomarch, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica, int transfEspecifica)
{
	L_ImageRGBUchar im;
	dibujaTransformaciones(im, calces1_o_caja0, maxImRef, maxTransfPorImagen, maxTransfTotales, imRefEspecifica, transfEspecifica);
#if __COMPAT_ATLIMAGE__
	return im.writeImageATL(nomarch)!=0;
#else
	return im.writeBMP(nomarch)!=0;
#endif
}

void L_RecPuntInt::imprimeReporte(FILE *fp)
{
	if (fp==NULL)
		fp=stdout;
	L_CountingTree contador;
	int i;
	fprintf(fp, "---------\n");
	fprintf(fp, "Reporte>>\n");
	fprintf(fp, "\n");
	fprintf(fp, "Calces iniciales totales: %ld\n", long(calL.size()));
	fprintf(fp, "Transformaciones totales encontradas: %ld\n", long(traL.size()));
	fprintf(fp, "\n");
	fprintf(fp, "Transformaciones por imagen de referencia:\n");
	traL.calcNumTransf(contador);
	for (i=0; i<nobj; i++)
		fprintf(fp, "   Transformaciones imagen referencia %d: %d\n", i, contador.votosDe(i));
	fprintf(fp, "\n");
	fprintf(fp, "Lista de todas las transformaciones:\n");
	traL.escribeResumen(fp);
	fprintf(fp, "--------\n");
}

void L_RecPuntInt::transformaImagenPruebaEnEntrenamiento()
{
	L_DescriptorLista desLAnt;
	L_ImageGrayDouble imAnt;
	desLPru.moveListTo(desLAnt);
	imAnt = *impru;
	reseteaPrueba();
	reseteaEntren();
	agregaImagenEntren_interno(imAnt, 0, 0.5, &desLAnt);
	finalizaEntren();
}

bool L_RecPuntInt::leeParametrosPropaga(const char *arch)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::leeParametrosPropaga");
	L_ParamLabelList l;
	FILE *fp;
	paramNoUsados.clear();
	fp=fopen(arch,"r");
	if (l.readFile(fp)==false)
	{
		L_POP_EXECUTING_FN("L_RecPuntInt::leeParametrosPropaga");
		return false;
	}
	paramRPI.updateValues_me_and_children(l);
	l.moveListTo(paramNoUsados);
	L_POP_EXECUTING_FN("L_RecPuntInt::leeParametrosPropaga");
	return true;
}

bool L_RecPuntInt::grabaParametrosPropaga(const char *arch)
{
	L_PUSH_EXECUTING_FN("L_RecPuntInt::grabaParametrosPropaga");
	L_ParamLabelList l;
	FILE *fp;
	paramRPI.buildListOfValues_me_and_children(l);
	paramNoUsados.copyListOn(l);
	fp=fopen(arch,"w");
	L_POP_EXECUTING_FN("L_RecPuntInt::grabaParametrosPropaga");
	return l.saveFile(fp);
}

void L_RecPuntInt::swap(L_RecPuntInt &other)
{
	other.L_RecPuntInt_subObjCopiable::swap(*this);
	paramRPI.swap(other.paramRPI);
	paramNoUsados.swap(other.paramNoUsados);
	corrPixBN.swap(other.corrPixBN);
	imsref.swap(other.imsref);
	calL.swap(other.calL);
	desLRef.swap(other.desLRef);
	desLPru.swap(other.desLPru);
	traL.swap(other.traL);
}



int L_RecPuntIntConjuntoTracks::procesaImagen(const L_ImageGrayDouble &im, const L_ImageGrayDouble &imc2, long frame)
{
	L_DescriptorNodo *desn;
	L_Calcec2Nodo *caln;
	L_DescriptorLista desLOri, desLc2;
	L_CalceLista calLOri, calLc2;
	L_Calcec2Lista calc2L;
	L_ImageGrayDouble imx;
	double tiempito = L_TIME(), dist; //, dt;
	int i, arrn;
	int nDes;
	bool acept, aceptc2;

	if (rec.gtc == NULL || rec.gtc->tipo != genTrAfinHough)
	{
		printf("L_RecPuntIntConjuntoTracks::procesaImagen() : object gtc del reconocedor mal definido\n");
		throw L_ArgException();
	}

	if (rec.imsref.size() == 0) // Primera llamada
	{
		ultId = 0;
		imx = im;
		rec.agregaImagenEntrenCopiaObjetoDe(imx);
		((L_GenTransfAfinHough *)rec.gtc)->hash.dxBin = 0.1;
		((L_GenTransfAfinHough *)rec.gtc)->hash.dyBin = 0.1;
		nDes = (int)rec.desLRef.size();
		if (rec.desLRef.size() == 0)
		{
			rec.reseteaEntren();
			return nDes;
		}
		rec.desLRef.sort(0, &L_Descriptor::qCmpIntens_nodo);
		for (desn = rec.desLRef.root; desn!=NULL; desn=desn->sig)
		{
			if ((int)arr.size() > nMax)
				break;
			arr.resize(arr.size()+1);
			arr[arr.size()-1].t.resize(0);
			arr[arr.size()-1].createFrom(desn->c, frame, ultId++,255,255,255);
		}
		rec.finalizaEntren();
		return nDes;
	}

	// "Envejecer" los tracks anteriores
	for (i=0; i<(int)arr.size(); i++)
	{
		// dt = tiempito - arr[i].t[arr[i].t.size()-1].ti;  // No se usa T_T

		if (frame - arr[i].t[ arr[i].t.size() - 1 ].frame > frameMax)
		{
			arr[i] = arr[arr.size()-1]; // Se sobrescribe el elemento actual
			arr.resize(arr.size()-1);
			i--; // Para descontar esta iteracion del for
		}
	}
	arrn = arr.size(); // arr.size() antes de comenzar a iterar, se modifica entre medio

	// Generar deteccion hacia im
	imx = im;
	rec.procesaImagenesPrueba(imx);
	if (rec.traL.size() == 0)
	{
		rec.transformaImagenPruebaEnEntrenamiento();
		return (int)rec.desLPru.size();
	}
	// Guardar deteccion de im
	rec.desLPru.moveListTo(desLOri);
	rec.traL.root->c->calL.moveListTo(calLOri);
	printf(" traL.size()=%ld: - - - - - - - - - - - - - - - - - ", long(rec.traL.size()));
	rec.reseteaPrueba();

	// Generar deteccion hacia im2
	if (imc2.data() != NULL)
	{
		imx = imc2;
		rec.procesaImagenesPrueba(imx);
		nDes = (int)rec.desLPru.size();
		// Guardar deteccion de im2
		rec.desLPru.moveListTo(desLc2);
		if (rec.traL.size() > 0)
			rec.traL.root->c->calL.moveListTo(calLc2);
		rec.reseteaPrueba();
	}

	printf("calLc2.size()=%ld\n", long(calLc2.size()));

	// Parear calces
	if (calLc2.size() > 0)
		calc2L.parear(calLOri, calLc2);
	else
		calc2L.agregar(calLOri);

	nDes = (int)desLOri.size();
	if ((int)desLc2.size() > nDes)
		nDes = (int)desLc2.size();

	// Procesar calces
	for (caln = calc2L.root; caln!=NULL; caln=caln->sig)
	{
		for (i=0; i<arrn; i++)
		{
			if (caln->c.dRef->id == arr[i].d.id && caln->c.dRef->tipoD == arr[i].d.tipoD &&
				caln->c.dRef->tipoP == arr[i].d.tipoP && arr[i].t[arr[i].t.size()-1].frame < frame)
			{
				if (caln->c.dPru != NULL)
					dist = (caln->c.dRef->x0 - caln->c.dPru->x0)*(caln->c.dRef->x0 - caln->c.dPru->x0) +
					(caln->c.dRef->y0 - caln->c.dPru->y0)*(caln->c.dRef->y0 - caln->c.dPru->y0);
				else
					dist = 1e100;
				dist = sqrt(dist);

				// He aqui la parte complicada
				// Se comparan los puntos respecto a los ultimos que se vieron
				acept = (dist < 30);

				if (caln->c.dPruc2 != NULL)
					dist = (caln->c.dRef->x0 - caln->c.dPruc2->x0)*(caln->c.dRef->x0 - caln->c.dPruc2->x0) +
					(caln->c.dRef->y0 - caln->c.dPruc2->y0)*(caln->c.dRef->y0 - caln->c.dPruc2->y0);
				else
					dist = 1e100;
				dist = sqrt(dist);
				aceptc2 = (dist < 30);

				switch (acept + 2*aceptc2)
				{
				case 0:
					break;
				case 1:
					const_cast<L_Descriptor *>(caln->c.dPru)->nobj = 10000;  // acept
					const_cast<L_Descriptor *>(caln->c.dPru)->id = caln->c.dRef->id;
					arr[i].addFrom(*caln->c.dPru, frame, tiempito);
					break;
				case 2:
					const_cast<L_Descriptor *>(caln->c.dPruc2)->nobj = 10000; // aceptc2
					const_cast<L_Descriptor *>(caln->c.dPruc2)->id = caln->c.dRef->id;
					arr[i].agrega_c2(*caln->c.dPruc2, frame, tiempito);
					break;
				case 3:
					const_cast<L_Descriptor *>(caln->c.dPru)->nobj = 10000; // acept && aceptc2
					const_cast<L_Descriptor *>(caln->c.dPru)->id = caln->c.dRef->id;
					const_cast<L_Descriptor *>(caln->c.dPruc2)->nobj = 10000;
					const_cast<L_Descriptor *>(caln->c.dPruc2)->id = caln->c.dRef->id;
					arr[i].addFrom(*caln->c.dPru, *caln->c.dPruc2, frame, tiempito);
					break;
				}
				break;
			}
		}
	}
	for (desn = desLOri.root; desn!=NULL; desn=desn->sig)
	{
		if (desn->c.nobj == 10000)
		{
			desn->c.nobj = 1;
			continue;
		}
		if ((int)arr.size() > nMax)
			break;
		arr.resize(arr.size()+1);
		arr[arr.size()-1].t.resize(0);
		arr[arr.size()-1].createFrom(desn->c, frame, ultId++,255,255,255);
	}

	// Esto puede hacer que los puntos "bailen", lo cual es algo molesto
	rec.desLPru.clear();
	for (i=0; i<(int)arr.size(); i++)
		rec.desLPru.push_back(arr[i].d); // La wea chanta pero es lo que hay
	rec.transformaImagenPruebaEnEntrenamiento();

	return nDes;
}










void L_RecPuntIntConjuntoTracks::pide_xy_n_frames(int nFrames, int &nTracks, L_Matrix &xy, long frame)
{
	int i, j, nel;
	bool conti;
	nTracks = 0;

	for (i=0; i<(int)arr.size(); i++)
	{
		if ((int)arr[i].t.size() < nFrames )
			continue;
		conti = false;
		nel = (int)arr[i].t.size();
		for (j=0; j<nFrames; j++)
		{
			if (arr[i].t[nel-nFrames+j].frame != frame - (nFrames-1) + j)
			{
				conti = true;
				break;
			}
		}
		if (conti)
			continue;
		nTracks++; // Este track sirve
	}
	if (nTracks==0)
	{
		xy.li = 0;
		xy.lj = 0;
		xy.destroy();
		return;
	}

	xy.reallocate(nTracks, 2*nFrames);
	nTracks = 0;

	for (i=0; i<(int)arr.size(); i++)
	{
		if ((int)arr[i].t.size() < nFrames )
			continue;
		conti = false;
		nel = arr[i].t.size();
		for (j=0; j<nFrames; j++)
		{
			if (arr[i].t[nel-nFrames+j].frame != frame - (nFrames-1) + j)
			{
				conti = true;
				break;
			}
		}
		if (conti)
			continue;
		for (j=0; j<nFrames; j++)
		{
			xy(nTracks,2*j+0)=arr[i].t[nel-nFrames+j].x;
			xy(nTracks,2*j+1)=arr[i].t[nel-nFrames+j].y;
		}
		nTracks++;
	}
}

void L_RecPuntIntConjuntoTracks::pide_uv_n_frames(int n, int &nTracks, const L_Camara *c, L_Matrix &uv, long frame)
{
	pide_xy_n_frames(n, nTracks, uv, frame);
	L_RayoCamara r;
	int i, j;
	for (i=0; i<nTracks; i++)
	{
		for (j=0; j<n; j++) // Esto se me habia ido, obvio
		{
			c->calcRayo(r, uv(i,2*j+0), uv(i,2*j+1));
			uv(i,2*j+0)=r.tanIzq;
			uv(i,2*j+1)=r.tanArr;
		}
	}
}


void L_RecPuntIntConjuntoTracks::pide_xy_n_frames_c2(int nFrames, int &nTracks, L_Matrix &xy, long frame)
{
	int i, j, nel;
	bool conti;
	nTracks = 0;

	for (i=0; i<(int)arr.size(); i++)
	{
		if ((int)arr[i].t.size() < nFrames )
			continue;
		conti = false;
		nel = arr[i].t.size();
		for (j=0; j<nFrames; j++)
		{
			if (arr[i].t[nel-nFrames+j].framec2 != frame - (nFrames-1) + j)
			{
				conti = true;
				break;
			}
		}
		if (conti)
			continue;
		nTracks++; // Este track sirve
	}
	if (nTracks==0)
	{
		xy.li = 0;
		xy.lj = 0;
		xy.destroy();
		return;
	}

	xy.reallocate(nTracks, 2*nFrames);
	nTracks = 0;

	for (i=0; i<(int)arr.size(); i++)
	{
		if ((int)arr[i].t.size() < nFrames )
			continue;
		conti = false;
		nel = arr[i].t.size();
		for (j=0; j<nFrames; j++)
		{
			if (arr[i].t[nel-nFrames+j].framec2 != frame - (nFrames-1) + j)
			{
				conti = true;
				break;
			}
		}
		if (conti)
			continue;
		for (j=0; j<nFrames; j++)
		{
			xy(nTracks,2*j+0)=arr[i].t[nel-nFrames+j].xc2;
			xy(nTracks,2*j+1)=arr[i].t[nel-nFrames+j].yc2;
		}
		nTracks++;
	}
}



//

L_RecPuntIntTrackPuntos::L_RecPuntIntTrackPuntos()
{
	nf = 5;
	ultTrack = 0;
	frame = 0;
	//
	sift.calces_maxRatio = 0.9; // Genrar bastantes calces, era 0.5 en el caso anterior
	//
	//tree.nCompMax = 64;
	//tree.usePercentage = false;
	//
	//hough.fusion.activo = false;  // Al parecer no influye en el tiempo
	hough.nMaxCeldasRev = 5;
	hough.mostrarTiemposInterno = false;
	//
	hough.hash.dxBin = 0.1;
	hough.hash.dyBin = 0.1;
	hough.hash.dangBin = 15*M_PI/180;
	hough.hash.doctBin = 0.5;
	hough.ransac.activo = true; // Hay que probar si funciona
	hough.ransac.dx = 0.6;
	hough.ransac.dy = 0.6;
	hough.topDown = false;
	hough.debFiltr.activo = false;
	hough.fusion.activo = false;
	//
	rec.corrPixBN.verif = false;
	rec.corrPixBN.debFiltr.activo = false;
	rec.RSL.apply = false; // Muy lnto
	rec.RSL.corrige = true;
}


void L_RecPuntIntTrackPuntos::main_procesar()
{
	double t1, t2;
	t1 = L_TIME();
	//L_CapturadorImagen capt("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\heman.avi");
	L_CapturadorImagen capt("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\KBRealSequence.avi");
	//L_CapturadorImagen capt("C:\\users\\ploncomi\\documents\\nti.avi");
	L_RecPuntIntTrackPuntos tracker;
	tracker.procesar(capt, 150);
	t2 = L_TIME();
	printf("Tiempo total: %f seg\n", t2-t1);
}

bool L_RecPuntIntTrackPuntos::procesar(L_CapturadorImagen &capt, int nMax)
{
	int i = 0;
	L_ImageGrayDouble imBN;
	while((nMax == -1 || i < nMax))
	{
		printf("Leyendo archivo... ");
		if (capt.capturarImagen() == false)
			break;
		printf("ok\n");
		imBN = capt.im;
		if (procesar_img(imBN) == false)
		{
			printf("Fallo en la iteracion %d", i);
			return false;
		}
		i++;
	}
	preparar_matriz();
	FILE *fp;
	fp = fopen("matches.txt", "w");
	if (fp == NULL)
	{
		printf("No se pudo crear el archivo matches.txt");
		return false;
	}
	fprintf(fp, "%d %d\n", ultTrack, frame); // ultTrack = num descr,  frame = num frames
	tracks.saveFile(fp);
	fclose(fp);
	return true;
}

bool L_RecPuntIntTrackPuntos::procesar_img(L_ImageGrayDouble &imNueva)
{
	L_MatrixBaseWidthStep<int> hip;
	std::vector<L_Descriptor *> desArr;
	L_DescriptorLista desLtmp;
	L_CalceNodo *cal;
	int i, k;
	// Redimensionar objetos
	imAnt.resize(nf);
	desAnt.resize(nf);
	// Formar reconocedor
	rec.fijaReconocedor(&sift, &tree, &hough);

	// Mover valores anteriores, [0] = actual, [nf-1] = el mas viejo

	imAnt[nf-1].destroy();
	desAnt[nf-1].clear();
	for (i=nf-1; i>=1; i--)
	{
		imAnt[i-1].swap(imAnt[i]);
		desAnt[i-1].moveListTo(desAnt[i]);
		desAnt[i].fijaPunteroImagen(&imAnt[i]);
	}

	imAnt[0] = imNueva;

	// Calcular descriptores de la imagen actual
	printf("Frame %d : calculando descriptores ", frame);
	sift.nobj = rec.nobj;
	sift.calcDescrPuntInt(imAnt[0]);
	printf("(%ld)\n", long(sift.desL.size()));

	printf("Frame %d : procesando calces\n", frame);
	// Enumerar descriptores de la imagen actual para poder referirse a ellos
	sift.desL.fijaIdCorrelativos();

	// Copiar descriptores actuales a los arreglos
	sift.desL.copyListOn(desAnt[0]);
	desAnt[0].copia_a_arreglo_punteros(desArr);

	//Crear matriz de hipotesis
	hip.reallocate_fc((int)desAnt[0].size(), nf);
	hip.setConstant(-1);

	// Generar kd-tree
	rec.agregaImagenesEntren_usa_desL(sift.desL, &imAnt[0], true);
	rec.finalizaEntren();

	// Generar correspondencias con imagenes anteriores, hip[0][:] no se usa
	for (k=1; k<nf; k++)
	{
		for (i=0; i<hip.li; i++)
			hip(i,k) = -1;
		if (imAnt[k].data() == NULL)
			continue;
		desAnt[k].copyListOn(desLtmp);
		rec.procesaImagenesPrueba_usa_desL(desLtmp, &imAnt[k]); // desLtmp es absorbido
		if (rec.traL.size() > 0)
		{
			for (cal = rec.traL.root->c->calL.root; cal!=NULL; cal=cal->sig) // dRef->id es el id del descriptor actual
				hip(cal->c.dRef->id,k) = cal->c.dPru->id; // dPru->id es el id del track
		}
		rec.reseteaPrueba();
	}
	//rec.ti.escribeResumen();
	rec.reseteaEntren();

	printf("Analizando hipotesis... ");
	// Analizar la matriz de hipotesis
	int hipCero, nErr;
	for (i=0; i<hip.li; i++)
	{
		hipCero = -1; // Numero de track actual
		nErr = 0;
		for (k=1; k<nf; k++)
		{
			if (hipCero == -1)
				hipCero = hip(i,k); // Podria ser zero de nuevo pero no afecta
			else
			{
				if (hipCero != hip(i,k) && hip(i,k)!=-1)
					nErr++;
			}
		}
		if (nErr*1.0/nf < 0.3 && hipCero != -1) // Aceptado
			desArr[i]->id = hipCero;
		else
			desArr[i]->id = ultTrack++;
	}

	if (tracks.li < (ultTrack+1) || tracks.lj < 3*(frame+1))
		tracks.conservativeResize(2*(ultTrack+1)+10, 2*3*(frame+1)+10);
	for (i=0; i<hip.li; i++)
	{
		// O son tracks antiguos, o son tracks nuevos
		tracks(desArr[i]->id,3*frame+0) = desArr[i]->x0;
		tracks(desArr[i]->id,3*frame+1) = desArr[i]->y0;
		tracks(desArr[i]->id,3*frame+2) = 1.0;
	}
	printf("ok\n");
	printf("Frame %d : ok\n", frame);
	frame++;
	return true;
}

void L_RecPuntIntTrackPuntos::preparar_matriz()
{
	int i, j;
	int xtr, ntr;
	tracks.conservativeResize(ultTrack, frame*3);

	// Eliminar tracks con un solo punto
	ntr = 0;
	for (i=0; i<ultTrack; i++)
	{
		xtr = 0;
		for (j=0; j<frame; j++)
		{
			if (tracks(i,3*j+0) != 0 || tracks(i,3*j+1) != 0 || tracks(i,3*j+2) != 0)
				xtr++;
		}
		if (xtr < 2)
			continue;
		for (j=0; j<frame; j++)
		{
			tracks(ntr,3*j+0) = tracks(i,3*j+0);
			tracks(ntr,3*j+1) = tracks(i,3*j+1);
			tracks(ntr,3*j+2) = tracks(i,3*j+2);
		}
		ntr++;
	}

	tracks.conservativeResize(ntr, 3*frame);
	ultTrack = ntr;

	// Convertir [0 0 0) en [-1 -1  1)
	for (i=0; i<ultTrack; i++)
	{
		for (j=0; j<frame; j++)
		{
			if (tracks(i,3*j+0) == 0 && tracks(i,3*j+1) == 0 && tracks(i,3*j+2) == 0)
			{
				tracks(i,3*j+0) = -1;
				tracks(i,3*j+1) = -1;
				tracks(i,3*j+2) = 1;
			}
		}
	}
	return;
}

void L_RecPuntIntTrackPuntos::resetear()
{
	int i;
	for (i=0; i<(int)imAnt.size(); i++)
		imAnt[i].destroy();
	for (i=0; i<(int)desAnt.size(); i++)
		desAnt[i].clear();
	tracks.destroy();
	frame = 0;
	ultTrack = 0;
}


void L_RecPuntIntTrackPuntos::genDibujosTracks()
{
	//L_CapturadorImagen capt("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\heman.avi");
	L_CapturadorImagen capt("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\KBRealSequence.avi");
	//L_CapturadorImagen capt("C:\\users\\ploncomi\\documents\\nti.avi");
	FILE *fp = fopen("matches.txt", "r");

	L_Matrix m;
	L_ImageGrayDouble im;
	L_ShapeArray lins;
	char buf[1000];
	int nf = 5;

	if (fgets(buf, 998, fp) == NULL)
		{printf("L_RecPuntIntTrackPuntos::genDibujosTracks() : formato de archio incorrecto\n"); return;}
	printf("Leyendo archivo %s\n", buf);
	m.readFile(fp);
	fclose(fp);

	// Dibujar tracks
	int t, d, k, k0;
	int x0, y0, x1, y1;
	
	for (t=0; t<m.lj/3; t++)
	{
		printf("Procesando imagen %d de %d\n", t, m.lj/3);
		lins.resize(0);
		if (capt.capturarImagen() == false)
			break;
		for (d=0; d<m.li; d++)
		{
			if (m(d,t*3+0) != -1 && m(d,t*3+1) != -1)
			{
				k0 = t-nf;
				if (k0 < 1)
					k0 = 1;
				for (k=k0; k<=t; k++)
				{
					x0 = (int)m(d, (k-1)*3+0);
					y0 = (int)m(d, (k-1)*3+1);
					x1 = (int)m(d, k*3+0);
					y1 = (int)m(d, k*3+1);
					if (x0 != -1 && y0 != -1 && x1 != -1 && y1 != -1)
						lins.drawLine(x0, y0, x1, y1, 255-(k-k0)*20, 255-(k-k0)*20, 255-(k-k0)*20);
				}
			}
		}
		capt.im.genDrawing(lins);
		sprintf(buf, "imsal%04d.png", t);
		capt.im.saveImage(buf);
		lins.resize(0);
	}
}











L_CalibradorCamaraRPI::L_CalibradorCamaraRPI()
{
	genProy.usar_4calces=true;
	genProy.radioError = 10;
	genProy.considerarCam = false;
	genProy.minCalces = 20;

	rpi.RSL.apply=false;
	rpi.RSL.corrige=false; // Va a sacar los poquitos calces que queden...

	lxMM = 188; // De puro flojo, es la imagen que uso...
	lyMM = 141;

	rpiEntrenado = false;

	rpi.fijaReconocedor(&genPuntos, &tree, &genProy);
}

bool L_CalibradorCamaraRPI::procesaImagenReferencia()
{
	if (rpiEntrenado == false)
	{
		if (imRef.data()==NULL || imRef.lx == 0 || imRef.ly==0)
		{
			printf("L_CalibradorCamara_unPlano::procesaNuevaImagen() : Intento de procesar imagen de referencia nula\n");
			return false;
		}
		rpi.reseteaEntren();
		rpi.agregaImagenEntrenCopiaObjetoDe(imRef);
		rpi.finalizaEntren();
		rpiEntrenado=true;
	}
	return true;
}

bool L_CalibradorCamaraRPI::procesaNuevaImagen(L_ImageGrayDouble &im)
{
	L_CalceNodo *cal;
	double u, v;
	if (im.data() == NULL || im.lx==0 || im.ly==0)
		return false;
	if (procesaImagenReferencia() == 0) // Por si no se ha procesado; verifica solo si es necesario hacerlo
		return false; // Imagen de referencia no definida
	if (genProy.considerarCam == false)
		camInter.corrigeImagenBilineal(im, imInter); // im se corrige para que funcione el reconocedor de puntos
	else
	{
		genProy.cam = camInter;
		imInter = im;
	}
	rpi.reseteaPrueba();
	rpi.procesaImagenesPrueba(imInter); // Se procesa la imagen enderezada
	if (rpi.traL.size()==0)
		return false;
	puntos.agregarCuadro();
	if (genProy.considerarCam == false)
	{
		for (cal = rpi.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
		{
			// Las coordenadas en pru estan corregidas
			camInter.distorsionaPixel(cal->c.dPru->x0, cal->c.dPru->y0, u, v); // Las coordenadas se distorsionan
			puntos.agregarPunto(cal->c.dRef->x0*lxMM/imRef.lx , cal->c.dRef->y0*lyMM/imRef.ly , 0.0 , u , v);
		}
	}
	else
	{
		//Las coordenadas en pru no estan corregidas
		for (cal = rpi.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
			puntos.agregarPunto(cal->c.dRef->x0*lxMM/imRef.lx , cal->c.dRef->y0*lyMM/imRef.ly , 0.0 , cal->c.dPru->x0 , cal->c.dPru->y0);
	}
	return true;
}

bool L_CalibradorCamaraRPI::grabarArchivoConfiguracion(const char *nomarch)
{
	FILE *fp;
	fp = fopen(nomarch,"w");
	if (fp==NULL)
		return false;
	fprintf(fp, "%s\n", nomImagen);
	fprintf(fp, "%lf %lf\n", lxMM, lyMM);
	fclose(fp);
	return true;
}

bool L_CalibradorCamaraRPI::leerArchivoConfiguracion(const char *nomarch)
{
	FILE *fp;
	size_t length;
	fp =fopen(nomarch, "r");
	if (fp==NULL)
		return false;
	if (fgets(nomImagen, 1000, fp) == NULL)
		return false;
	if (fscanf(fp, "%lf%lf", &lxMM, &lyMM) < 2)
		return false;
	length=strlen(nomImagen);
	if (nomImagen != 0)
	{
		while (nomImagen[length-1] == '\n')
			nomImagen[--length]=0;
		while (nomImagen[length-1] == ' ')
			nomImagen[--length]=0;
		if (imRef.readImage(nomImagen) == 0)
			return false;
	}
	return true;
}

bool L_CalibradorCamaraRPI::calibrar_Tsai()
{
	FILE *fp;
	double arr[20];
	double arrTmp[20];
	int sum=0;
	int i, j;

	for (j=0; j<20; j++)
	{
		arr[j] = 0;
		arrTmp[j]=0;
	}

	for (i=0; i<(int)puntos.size(); i++)
	{
		printf("Enviando %d de %ld\n", i, long(puntos.size()));
		puntos[i].guardar_TsaiCM("_temp1715.txt");
		if( system(L_DEL_SYSTEM_CMD " _temp1716.txt") ) {}
		if( system("TsaiCM ccal _temp1715.txt > _temp1716.txt") != 0)
		{
			printf("Error en TsaiCM en punto %d de %ld...\n", i, long(puntos.size()));
			continue;
		}
		fp = fopen("_temp1716.txt", "r");
		if (fp==NULL)
			continue;
		for (j=0; j<20; j++)
		{
			if (fscanf(fp, "%lf", &arrTmp[j]) < 1)
			{
				printf("Problema de lectura en _tem1716.txt linea %d\n", j);
			}
		}
		fclose(fp);
		if (arrTmp[19] > 30) // Error de ajuste muy grande
			continue;
		for (j=0; j<20; j++)
			arr[j]+=arrTmp[j];
		sum++;
	}

	for (j=0; j<20; j++)
		arr[j]/=sum;

	//    0          1          2           3            4            5          6        7        8
	// Ncx[sel], Nfx[pix], dx[mm/sel], dy[mm/sel], dpx[mm/pix], dpy[mm/pix], Cx[pix], Cy[pix], sx[df/df]
	//    9         10         11
	// f[mm/tan], k1[1/mm^2], k2[1/mm^4]
	double dfX, dfY;
	dfX = arr[8]/arr[4] * arr[9];  // sx/d'x * f
	dfY = 1.0   /arr[5] * arr[9]; //  1/d'y * f
	double f2 = arr[9]*arr[9];
	double k1 = arr[10]*f2; // k1 * f^2
	double k2 = 0;
	//camInter.fijaDistorsionK1K2(k1, k2);
	camInter.fijaParametros(camInter.resXv(), camInter.resYv(), camInter.cenYv(), camInter.cenYv(), k1, k2, 0, 0, dfX, dfY, 0);

	return true;
}



L_TemplateDescr::L_TemplateDescr() : gUnion(4), gMi(NULL,NULL,NULL), paramTempl("Templ", 4)
{
	compr = false;

	gUnion.agregaReferencia(&gMi);
	gUnion.agregaReferencia(&gLo);
	gUnion.agregaReferencia(&gHa);
	gUnion.agregaReferencia(&gHa1);
	
	paramTempl.addChildren(gMi.pideParams());
	paramTempl.addChildren(gLo.pideParams());
	paramTempl.addChildren(gHa.pideParams());
	paramTempl.addChildren(gHa1.pideParams());

	gMi.agregaEspejo = true;
}

void L_TemplateDescr::fijaGeneradores(bool usaMinuc, bool usaLSift, bool usaPHarris, bool usaHarris)
{
	gUnion.activo[0]=usaMinuc;
	gUnion.activo[1]=usaLSift;
	gUnion.activo[2]=usaPHarris;
	gUnion.activo[3]=usaHarris;
}

bool L_TemplateDescr::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	bool ret;
	ret = gUnion.calcDescrPuntInt(im, s0);
	gUnion.desL.moveListTo(desL);
	return ret;
}


bool L_TemplateDescr::procesaEntrada(const char *name)
{
	double t1;
	t1 = L_TIME();
	if (im.readImage(name)==0) // Imagefile
	{
		printf("No se puede abrir imagefile en ENROLL\n");
		return false;
	}
	if (calcDescrPuntInt(im, 0.5)==false)
	{
		printf("Problemas al calcular descriptores locales... :(\n");
		return false;
	}
	strcpy(nomIm, name);
	//int largoNomArch = (int)strlen(name);
	L_FileName d(name);
	d.ext = "txt";
	d.name = d.name + "_min";  //  "xxx.bmp"  ->  "xxx_min.bmp" , "xxx_min.txt"
	gMi.fijaArchivoMinucias(d.c_str());
	gMi.buscaCamposOrientacionFrecuencia();
	dt = L_TIME() - t1;
	return true;
}

bool L_TemplateDescr::leeTemplate(const char *name)
{
	FILE *fp;
	L_ImageGrayUchar imBN;
	int largoNomArch;

	fp=fopen(name,"rb"); // templatefile
	if (fp==NULL)
	{
		printf("No se puede abrir archivo template %s\n", name);
		return false;
	}
	// fp: desL, imBN, dt, largoNomArch, argv[1], camposOrientacionFrecuencia
	
	desL.clear(); // Por si acaso
	compr = (fgetc(fp) != 0);
	if (compr)
	{
		std::vector<char> bufM, bufC;
		int enee;
		if (fread(&enee, sizeof(enee), 1, fp) != 1)
			{fclose(fp); return false;}
		bufC.resize(enee);
		if (fread(&(bufC[0]), sizeof(char), bufC.size(), fp) != bufC.size())
			{fclose(fp); return false;}
		L_HuffmanEncoder::decodeAll(bufC, bufM);
		// gUnion se evaluate recursivamente para encontrar el refGen adecuado
		desL.readListMem(bufM, 0, &im, &gUnion); // Se le pasa la referencia aunque la imagen no se haya leido aun
	}
	else
		desL.readListBin(fp, &im, &gUnion); // Se uso imPru como puntero-imagen para los descriptores del template. Cuidado con esta linea

	imBN.leeComprGris(fp);
	im = imBN;
	int nleidos = 0;
	nleidos += (int)fread(&dt,sizeof(dt),1,fp);
	nleidos += (int)fread(&largoNomArch, sizeof(largoNomArch), 1, fp);
	nleidos += (int)fread(nomIm, 1, L_MIN(largoNomArch+1, 999), fp);
	if (nleidos < 1 + 1 + largoNomArch+1)
		return false;
	gMi.leeCamposOrientacionFrecuencia(fp);
	fclose(fp); // cierra templatefile
	return true;
}

bool L_TemplateDescr::grabaTemplate(const char *name)
{
	int largoNomArch;
	L_ImageGrayUchar imBN;
	FILE *fp;

	fp=fopen(name,"wb"); // Template
	if (fp==NULL)
	{
		printf("No se puede generar archivo template en ENROLL\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
		return 1;
	}
	// fp: desL, imBN, dt, largoNomArch, argv[1], camposOrientacionFrecuencia

	// Primero: ver si se va a comprimir el resultado o no
	fputc(compr ? 1 : 0, fp);
	if (compr)
	{
		std::vector<char> bufM, bufC;
		int enee;
		desL.writeListMem(bufM);
		//L_HuffmanEncoder::encodeDelta(bufM, bufC);
		//L_HuffmanEncoder::encodeRLE0(bufC, bufM);
		L_HuffmanEncoder::encodeAll(bufM, bufC);
		fwrite(&enee, sizeof(enee), 1, fp);
		bufC.resize(enee);
		fwrite(&(bufC[0]), sizeof(char), bufC.size(), fp);
	}
	else
		desL.writeListBin(fp);

	imBN = im;
	imBN.grabaComprGris(fp, NULL, compr ? -1 : 0);

	largoNomArch = (int)strlen(nomIm);
	fwrite(&dt,sizeof(dt),1,fp);
	fwrite(&largoNomArch, sizeof(largoNomArch), 1, fp);
	fwrite(nomIm, 1, largoNomArch+1, fp);

	gMi.grabaCamposOrientacionFrecuencia(fp, compr);

	fclose(fp); // cierre Template
	return true;
}


/*
int L_ProgramasRecPuntInt::main_MuestraPuntos(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_MuestraPuntos");
	L_ImageGrayDouble imx;
	L_ImageRGBUchar imRGB;
#ifdef __L_LOWE_H_
	L_LoweGenDescr gpi;
#else
	L_PDoG_SIFT gpi;
#endif
	L_ShapeArray monos;
	if (argc<3 || argc>4 || argv[1][0]=='/')
	{
		printf("\n\n  muestraPuntos imagen.bmp output.bmp porcLum\n\n");
		printf("porcLum: numero entre 1 y 100 que indica la luminosidad mean de la output\n\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_MuestraPuntos");
		return 1;
	}
	if (imx.readBMP(argv[1])==0)
	{
		printf("No se encuentra el archivo\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_MuestraPuntos");
		return 1;
	}
	double porc = 100.0;
	if (argc == 4)
		sscanf(argv[3], "%lf", &porc);
	L_Lowe_Generator::DoubleImSize=false;
	gpi.calcDescrPuntInt(imx, 0.5);
	gpi.desL.dibujaFlechas(monos,255,255,255,2);
	imx.multiply_to_each_element(porc/100.0);
	imRGB=imx;
	imRGB.genDrawing(monos);
	imRGB.writeBMP(argv[2]);
	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_MuestraPuntos");
	return 0;
}
*/

int L_ProgramasRecPuntInt::main_GrafPuntInt(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
	L_ImageGrayDouble imPru;
	L_ImageGrayDouble imRef;
	L_Array<const L_ImageGrayDouble *> imarr;
	L_ImageRGBUchar imSal;
	L_ShapeArray lins;
	char arch[100];
	char tipo[100];
	char lin[1000];
	FILE *fp;
	int ret1=1, ret2=1;
#ifdef __ATLIMAGE_H__
	L_ImageRGBUchar imRGB;
#endif

	if (argc!=6)
	{
		printf("\n\n  GrafPuntInt tipo archCal.txt imRef imPru imSal\n\n");
		printf("    Tipo: nTrans, mejoresCalces, mejorTransf\n\n");
#ifdef __ATLIMAGE_H__
		printf("Las imagenes pueden ser .bmp, .jpg u otras extensiones\n\n");
#else
		printf("Las imagenes solo pueden ser .bmp\n\n");
#endif
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
		return 1;
	}
	sprintf(tipo, ">>>>:%s\n", argv[1]);
	sprintf(arch, ">>Pru: %s\n", argv[4]);
	fp=fopen(argv[2],"r");
	if (fp==NULL)
	{
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
		return 1;
	}
	if (*argv[1]!='n') // Hay que dibujar
	{
#ifdef __ATLIMAGE_H__
		ret1=imRGB.readImageATL(argv[3]);

		if (ret1!=0)
			imRef=imRGB;
		ret2=imRGB.readImageATL(argv[4]);
		if (ret2!=0)
			imPru=imRGB;

#else
		ret1=imRef.readBMP(argv[3]);
		ret2=imPru.readBMP(argv[4]);
#endif
		if (ret1==0 || ret2==0)
		{
			L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
			printf("Error: Imagen '%s' no existente\n", (ret1==0)?argv[3]:argv[4]);
			return 1;
		}
	}
	while (fgets(lin,999,fp)!=NULL)
		if (strcmp(lin,arch)==0)
			break;
	while (fgets(lin,999,fp)!=NULL)
		if (strcmp(lin,tipo)==0)
			break;
	if (feof(fp))
	{
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
		return 1;
	}
	if (*argv[1]!='n')
	{
		lins.read(fp);
		imarr.resize(2);
		imarr[0] = &imPru;
		imarr[1] = &imRef;
		imSal.genDrawingMatches(imarr, lins);
		imSal.writeBMP(argv[5]);
		fclose(fp);
	}
	else
	{
		if (fscanf(fp, "%d", &ret1) != 1)
			printf("Formato de archivo incorrecto\n");
		printf("%d\n",ret1);
		fclose(fp);
	}
	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_GrafPuntInt");
	return 0;
}

int L_ProgramasRecPuntInt::main_RecPuntInt(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");

	L_RecPuntInt rec; // reconocedor generico
	//L_SDoG_SIFT gpi; // generador de descriptores - puntos de interes
	//L_SHarris_SIFT gpi;
#ifdef __L_LOWE_H_
	L_LoweGenDescr gpi;
#else
	L_PDoG_SIFT gpi;
#endif
	L_KdTreeBBF busqArbol;  // generador de calces
	//L_BusqPermPiv busqPiv;  // generador de calces
	L_GenTransfAfinHough gtc;   // generador de transformaciones afines 2D a partir de calces

	L_ImageGrayDouble imRef;
	L_ImageGrayDouble imPru;
	L_ShapeArray lins;
	int i, ret;
	FILE *fp;
#ifdef __ATLIMAGE_H__
	L_ImageRGBUchar imRGB;
#endif

	if (argc<5)
	{
		printf("\n\n  RecPuntInt archParam.txt archCal.txt imRef imPru1 ... imPruN\n\n");
#ifdef __ATLIMAGE_H__
		printf("Las imagenes pueden ser .bmp, .jpg u otras extensiones\n\n");
#else
		printf("Las imagenes solo pueden ser .bmp\n\n");
#endif
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");
		return 1;
	}
	rec.fijaReconocedor(&gpi, &busqArbol, &gtc); // El reconocedor va a usar estos objetos cuando lo necesite
	rec.leeParametrosPropaga(argv[1]);
	rec.grabaParametrosPropaga(argv[1]); // Si no hay arch de param, lo createFrom
	fp=fopen(argv[2],"w");

	if (fp==NULL)
	{
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");
		printf("No se puede crear el archivo de calces '%s'\n", argv[2]);
		return 1;
	}

#ifdef __ATLIMAGE_H__
	ret=imRGB.readImageATL(argv[3]);
	if (ret!=0)
		imRef=imRGB;
#else
	ret=imRef.readBMP(argv[3]);
#endif
	if (ret==0)
	{
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");
		printf("No se encuentra la imagen de referencia '%s'\n", argv[3]);
		return 1;
	}

	rec.agregaImagenEntrenCopiaObjetoDe(imRef);
	rec.finalizaEntren();

	for (i=0; i<=argc-5; i++)
	{
		fprintf(fp, ">>Ref: %s\n", argv[3]);
		fprintf(fp, ">>Pru: %s\n", argv[i+4]);
#ifdef __ATLIMAGE_H__
		ret=imRGB.readImageATL(argv[i+4]);
		if (ret!=0)
			imPru=imRGB;
#else
		ret=imPru.readBMP(argv[i+4]);
#endif
		if (ret==0)
		{
			L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");
			printf("No se encuentra la imagen de test '%s'\n", argv[i+4]);
			return 1;
		}
		rec.procesaImagenesPrueba(imPru, 0.5);
		fprintf(fp,">>>>:nTrans\n");
		fprintf(fp,"%ld\n", long(rec.traL.size()));
		fprintf(fp,">>>>:mejoresCalces\n");
		if (rec.traL.size()>0)
			rec.traL.root->c->calL.genLineas(lins);
		lins.write(fp);
		lins.resize(0);
		fprintf(fp,">>>>:mejorTransf\n");
		rec.traL.dibujaTransformaciones(lins, 0, 1, 1, 1);
		lins.write(fp);
		lins.resize(0);
		rec.traL.destroyList(&rec.cal_pool);
		printf("Imagen de test %d de %d (%d%%)\n", i+1, argc-4, (int)(i*100.0/(argc-4)));
		rec.ti.escribeResumen(stdout);
		printf("\n");
	}
	fclose(fp);

	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_RecPuntInt");
	return 0;
}

int L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
	if (argc!=6)
	{
		printf("\n\n  GenCalcesTransf num archParam.txt archCal.txt imRef.bmp imPru.bmp\n\n");
		printf("num: 0 para sDoG + SIFT,  1 para sHarris + SIFT");
		#ifdef __L_LOWE_H_
			printf(",  2 para Lowe original");
		#endif
		printf("\n\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
		return 1;
	}

	L_RecPuntInt rec; // reconocedor generico
	L_GenDescrPuntInt *gpi;
	L_SDoG_SIFT gpi0;
	L_SHarris_SIFT gpi1;

	#ifdef __L_LOWE_H_
		L_LoweGenDescr gpi2;
	#endif
	L_GenTransfAfinHough gtc;   // generador de transformaciones afines 2D a partir de calces
	L_KdTreeBBF busqArbol;  // generador de calces
	L_ImageGrayDouble imRef;
	L_ImageGrayDouble imPru;
	FILE *fp;

	if (argv[1][0]=='0')
		gpi=&gpi0;
	else if (argv[1][0]=='1')
		gpi=&gpi1;
	#ifdef __L_LOWE_H_
		else if (argv[1][0]=='2')
			gpi=&gpi2;
	#endif
	else
	{
		printf("Error en argumento num\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
		return 1;
	}

	fp=fopen(argv[3],"w");
	if (fp==NULL)
	{
		printf("Error en argumento archCal.txt\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
		return 1;
	}

	rec.fijaReconocedor(gpi, &busqArbol, &gtc);
	rec.leeParametrosPropaga(argv[2]);
	rec.grabaParametrosPropaga(argv[2]);
	if ( imRef.readBMP(argv[4]) == 0 )
	{
		printf("Error en argumento imRef\n");
		fclose(fp);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
		return 1;
	}
	rec.agregaImagenEntrenCopiaObjetoDe(imRef);
	rec.finalizaEntren();
	if ( imPru.readBMP(argv[5]) == 0 )
	{
		printf("Error en argumento imPru\n");
		fclose(fp);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
		return 1;
	}
	rec.procesaImagenesPrueba(imPru);
	if (rec.traL.size()>0)
	{
		L_CalceNodo *cal;
		for (cal=rec.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
			fprintf(fp, "%.1f %.1f %.1f %.1f\n", cal->c.dRef->x0, cal->c.dRef->y0, cal->c.dPru->x0, cal->c.dRef->y0);
	}
	fclose(fp);
	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_ProcesaImagenesGenTxtCalces");
	return 0;
}

int L_ProgramasRecPuntInt::analizaResultadosUch(const char *metodo, int n1, int n2, int n3, int n4, int width)
{
	int i, j;
	FILE *fp;
	long nParesTot=0;
	long nParesIg=0;
	long nCorr=0;
	long nIncorr=0;
	long nTransfCorr=0;
	long nTransfIncorr=0;
	int numero;
	char nomarch[100];
	char lin1[100], lin2[100];
	L_String sNum1, sNum2;

	printf("Analizando imagenes de %d:%d contra %d:%d, metodo %s\n", n1, n2, n3, n4, metodo);

	for (i=n1; i<=n2; i++)
	{
		for (j=n3; j<=n4; j++)
		{
			sNum1=L_String(i, width);
			sNum2=L_String( (j-1)/10 , width);
			sprintf(nomarch, "calN_%s_%s_%s.txt", metodo, sNum1.c_str(), sNum2.c_str());
			fp=fopen(nomarch,"r");
			if (fp==NULL)
				continue;
			sNum2=L_String(j, width);

			sprintf(lin1, ">>Pru: uch%sb.bmp\n",sNum2.c_str());
			while(fgets(lin2, 99, fp)!=NULL && strcmp(lin1, lin2)!=0)
				;
			sprintf(lin1, ">>>>:nTrans\n");
			while(fgets(lin2, 99, fp)!=NULL && strcmp(lin1, lin2)!=0)
				;
			if (fgets(lin2, 99, fp)==NULL)
				continue;
			if (sscanf(lin2, "%d", &numero)==1)
			{
				nParesTot++;
				if (i==j)
				{
					nParesIg++;
					if (numero>0)
					{
						nCorr++;
						nTransfCorr+=numero;
					}
				}
				else
				{
					if (numero>0)
					{
						nIncorr++;
						nTransfIncorr+=numero;
					}
				}
			}
			fclose(fp);
		}
	}
	printf("Analizados: %ld de %ld = %.2f%%\n", nParesTot, (n2-n1+1L)*(n4-n3+1L), nParesTot*100.0/((n2-n1+1L)*(n4-n3+1L)));
	printf("Existe transf corr: %.2f%% ; nTransfProm=%.2f\n", nCorr*100.0/nParesIg, nTransfCorr*1.0/nCorr);
	printf("Existe transf inc:  %.2f%% ; nTransfInc =%.2f\n\n", nIncorr*100.0/(nParesTot-nParesIg), nTransfIncorr*1.0/nIncorr);
	return 0;
}

int L_ProgramasRecPuntInt::main_VerificaResultados(int argc, char *argv[])
{
	if (argc!=6)
	{
		printf(
			"\n    VerifRes imRef.bmp imPru.bmp archDet.txt archReal.txt archSal.txt\n\n"
			"  Verifica que la transformacion principal detectada (en archDet.txt) sea compatible\n"
			"con la indicada por la base de datos (en archReal.txt).\n"
			"  Se addFrom un 0 o un 1 a archSal.txt para indicar si la deteccion fue exitosa\n\n"
			);
		return 1;
	}
	L_FileName dna;
	char *nomImRef=argv[1];
	char *nomImPru=argv[2];
	char *nomTxtDet=argv[3];
	char *nomTxtReal=argv[4];
	char *nomTxtSal=argv[5];
	FILE *fp;
	char lin[200], nomtmp[100], *ptrS;
	char nomImPru_[100], nomImRef_[100];
	L_GenTransfAfinHough gtah;
	L_TransfAfinPI2D trDet;
	L_TransfAfinPI2D trReal;
	double xPru, yPru, xRef, yRef;
	long nt1=0, nt2=0, nt3=0;
	L_DescriptorLista d1, d2;
	bool aceptar=true;
	bool res=false;
	L_ImageGrayDouble im;

	im.lx=640;
	im.ly=480;

	gtah.fusion.activo=true;
	gtah.fusion.superp=0.85;

	fp=fopen(nomTxtDet,"r");
	if (fp==NULL)
	{
		printf("Archivo de transf detectada \"%s\" no encontrado\n",nomTxtDet);
		return 1;
	}
	dna.asign(nomImPru);
	strcpy(nomImPru_,dna.name.data());
	strcat(nomImPru_,".");
	strcat(nomImPru_,dna.ext.data());
	while( (ptrS=fgets(lin,200,fp)) != NULL )
	{
		if ( strlen(lin)>6 && strncmp(lin, ">>Pru:",6)==0)
		{
			sscanf(lin,">>Pru: %s", nomtmp);
			if (strcmp(nomtmp,nomImPru_)!=0)
				continue;
			break; // Se logro encontrar la transformacion detectada
		}
	}
	if (ptrS==NULL)
	{
		printf(
			"Info de la imagen:\n"
			"  Pru: \"%s\" \n"
			"    no encontrada en %s\n", nomImPru, nomTxtDet);
		fclose(fp);
		return 1;
	}
	// Leer la info de la transformacion detectada
	while( (ptrS=fgets(lin,200,fp)) != NULL )
	{
		if ( strlen(lin)>6 && strncmp(lin, ">>Pru:",6)==0)
		{
			ptrS=NULL;
			break;
		}
		if ( strcmp(lin, ">>>>:mejoresCalces\n")!=0)
			continue;
		while (fscanf(fp, "%lf %lf %*d %lf %lf %*d %*d %*d %*d %*d", &xRef, &yRef, &xPru, &yPru)==4)
			trDet.calL.insertaNodoAlInicio_(d1,d2,(int)xRef, (int)yRef, (int)xPru, (int)yPru, &im, &im);
		break;
	}
	if (ptrS==NULL)
	{
		printf("Info de %s en %s no incluye mejor transformacion\n",nomImPru,nomTxtDet);
		return 1;
	}
	if (trDet.calcTransfMinCuad()==false) // No hay mejor transformacion o mejores calces
		aceptar=false;
	fclose(fp);
	trDet.calL.copyListOn(trReal.calL); // wea penosa que hay que hacer

	fp=fopen(nomTxtReal,"r");
	if (fp==NULL)
	{
		printf("Archivo de datos reales \"%s\" no encontrado\n",nomTxtReal);
		return 1;
	}
	dna.asign(nomImRef);
	strcpy(nomImRef_,dna.name.data());
	strcat(nomImRef_,".");
	strcat(nomImRef_,dna.ext.data());
	bool enArchivo=false;
	while( (ptrS=fgets(lin,200,fp)) != NULL )
	{
		if ( strlen(lin)>6 && strncmp(lin, ">>Ref:",6)==0)
		{
			sscanf(lin,">>Ref: %s", nomtmp);
			if (strcmp(nomtmp,nomImRef_)!=0)
				continue;
			if ( (ptrS=fgets(lin,200,fp)) == NULL )
			{
				printf("Error: linea >>Ref seguida por final de archivo en %s\n",nomTxtDet);
				fclose(fp);
				return 1;
			}
			if ( strlen(lin)<6 || strncmp(lin, ">>Pru:\n",6)!=0 )
			{
				printf("Error: linea >>Ref no seguida por linea >>Pru en %s\n",nomTxtDet);
				fclose(fp);
				return 1;
			}
			sscanf(lin,">>Pru: %s", nomtmp);
			if (strcmp(nomtmp,nomImPru)!=0)
			{
				enArchivo=false;
				continue;
			}
			enArchivo=true;
			continue; // Se logro encontrar la transformacion detectada
		}
		else if (enArchivo)
		{
			if ( strcmp(lin, ">>>>:TransfAfin\n")!=0)
				continue;
			ptrS=fgets(lin, 200, fp);
			if (ptrS==NULL || sscanf(lin, "%lf %lf %lf %lf %lf %lf", &trReal.m11, &trReal.m12, &trReal.tx, &trReal.m21, &trReal.m22, &trReal.ty) != 6)
			{
				printf("Parametros de la mejor transformacion incorrectos en %s\n", nomTxtDet);
				fclose(fp);
				return 1;
			}
			if (aceptar)
				res=gtah.verifFusion(trDet,trReal);
			break;
		}
	}
	fclose(fp);
	fp=fopen(nomTxtSal,"r");
	if (fp!=NULL)
	{
		if (fscanf(fp,"%ld%ld%ld",&nt1,&nt2,&nt3) != 3)
			printf("formato de archivo incorrecto\n");
		fclose(fp);
	}
	if (res)
		nt1++;
	else if (aceptar)
		nt2++;
	else
		nt3++;
	fp=fopen(nomTxtSal,"w");
	fprintf(fp,"%ld\n%ld\n%ld\n",nt1,nt2,nt3);
	fclose(fp);
	return 0;
}

int L_ProgramasRecPuntInt::main_CuentaDetecciones(int argc, char *argv[])
{
	if (argc!=3)
	{
		printf(
			"\n    CuentaDet archDet.txt archSal.txt\n\n"
			"  Cuenta las transformaciones principales detectadas en archDet.txt\n"
			"  Se addFrom un 0 o un 1 a archSal.txt para indicar si la deteccion fue exitosa\n\n"
			);
		return 1;
	}
	FILE *fp;
	char lin[200];
	int nt1=0, nt2=0;
	fp=fopen(argv[2],"r");
	if (fp!=NULL)
	{
		if (fscanf(fp,"%d%d",&nt1,&nt2) != 2)
		{
			printf("Formato de archivo incorrecto\n");
			fclose(fp);
			return 1;
		}	
		fclose(fp);
	}
	fp=fopen(argv[1],"r");
	if (fp==NULL)
	{
		printf("Error: Archivo %s no existe\n",argv[1]);
		return 1;
	}
	while(fgets(lin, 199, fp)!=NULL)
	{
		if (strcmp(lin,">>>>:mejoresCalces\n")==0)
		{
			if (fgets(lin, 199, fp)==NULL)
			{
				printf("Error en estructura de %s\n",argv[1]);
				fclose(fp);
				return 1;
			}
			if (strcmp(lin,"---\n")==0)
				nt2++; // No detectado
			else
				nt1++;
		}
	}
	fclose(fp);
	fp=fopen(argv[2],"w");
	if (fp==NULL)
	{
		printf("No se puede escribir en archivo %s\n",argv[2]);
		return 1;
	}
	fprintf(fp,"%d\n%d\n",nt1,nt2);
	fclose(fp);
	return 0;
}

////

int L_ProgramasRecPuntInt::main_FVC_enroll(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");

	double t1, t2, dt;

	t1=L_TIME();

	if (argc!=5)
	{
		printf(
			"\n"// 0      1           2          3           4
			"  ENROLL imagefile templatefile configfile outputfile\n"
			"    imagefile: una imagen de referencia\n"
			"    templatefile: archivo de output (datos de la imagen de referencia)\n"
			"    configfile: archivo.txt (parametros)\n"
			"    outputfile: archivo que indica la eficiencia del proceso\n"
			"Al usar COMPR=true se comprime (un poquito) el archivo template\n"
			"Al usar usaMINEXT=true, se pueden agregar los archivos:\n"
			"    imagefile_min.txt    (minucias)\n"
			"    imagefile_min_OF.bmp (campos de orientacion y frecuencia en canales R, G)\n"
			"Esos archivos se agregan al template para poder extraer caracteristicas extra\n"
			);
#ifndef __ATLIMAGE_H__
		printf("  --- Esta version read solo archivos BMP\n");
#else
		printf("  --- Esta version read varios formatos de imagen\n");
#endif
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
		return 1;
	}

	L_ImageRGBUchar imRGB;
	L_ImageGrayDouble imx;
	L_ImageGrayUchar imBN;

	int largoNomArch = (int)strlen(argv[1]);
	L_FileName d(argv[1]);
	d.ext = "txt";
	d.name = d.name + "_min";  //  "xxx.bmp"  ->  "xxx_min.bmp" , "xxx_min.txt"
	// Variables usadas para el reconocedor
	L_GenDescrPuntInt *gpi;
	L_MinuciasExterno gMi(d.c_str());
#ifdef __L_LOWE_H_
	L_LoweGenDescr gLo;
#else
	L_PDoG_SIFT gLo;
#endif
	L_PHarris_SIFT gHa;
	L_Harris_SIFT gHa1;
	L_GenDescrArr gUnion(4);
	L_ParamLabelList pli;
	L_ParamManagerLocal parLoc("FVC",5);
	// Parametros para controlar los generadores de descriptores usados
	bool usaSIFT=true;
	bool usaMINEXT=false;
	bool usaHARRIS=false;
	bool usaHARRIS1=false;

	bool COMPR = false;
	FILE *fp;

	// Configurar detectores de puntos
	gMi.agregaEspejo=true;
	gUnion.agregaReferencia(&gMi);
	gUnion.agregaReferencia(&gLo);
	gUnion.agregaReferencia(&gHa);
	gUnion.agregaReferencia(&gHa1);
	gpi=&gUnion; // Usar todos a la vez
	
	// Configurar parametros
	parLoc.addFrom("usaMINEXT",&usaMINEXT);
	parLoc.addFrom("usaSIFT",&usaSIFT);
	parLoc.addFrom("usaHARRIS",&usaHARRIS);
	parLoc.addFrom("usaHARRIS1",&usaHARRIS1);
	parLoc.addFrom("COMPR",&COMPR);
	parLoc.addChildren(gUnion.pideParams());
	pli.readFile(fopen(argv[3],"r")); // Configfile
	parLoc.updateValues_me_and_children(pli); // pli contiene ahora los parametros que sobraron
	parLoc.buildListOfValues_me_and_children(pli);
	pli.saveFile(fopen(argv[3],"w"));
	gUnion.activo[0]=usaMINEXT;
	gUnion.activo[1]=usaSIFT;
	gUnion.activo[2]=usaHARRIS;
	gUnion.activo[3]=usaHARRIS1;

	// Lectura de la imagen
	if (imRGB.readImage(argv[1])==0) // Imagefile
	{
		printf("No se puede abrir imagefile en ENROLL\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
		return 1;
	}

	imBN=imRGB;
	imx=imBN;

	// Calculo de puntos de interes
	if (gpi->calcDescrPuntInt(imx, 0.5)==false)
	{
		printf("Problemas al calcular descriptores locales... :(\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
		return 1;
	}

	fp=fopen(argv[2],"wb"); // Template
	if (fp==NULL)
	{
		printf("No se puede generar archivo template en ENROLL\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
		return 1;
	}
	// fp: desL, imBN, dt, largoNomArch, argv[1], camposOrientacionFrecuencia

	// Primero: ver si se va a comprimir el resultado o no
	fputc(COMPR ? 1 : 0, fp);
	if (COMPR)
	{
		std::vector<char> bufM, bufC;
		int enee;
		gpi->desL.writeListMem(bufM);
		//L_HuffmanEncoder::encodeDelta(bufM, bufC);
		//L_HuffmanEncoder::encodeRLE0(bufC, bufM);
		L_HuffmanEncoder::encodeAll(bufM, bufC);
		fwrite(&enee, sizeof(enee), 1, fp);
		bufC.resize(enee);
		fwrite(&(bufC[0]), sizeof(char), bufC.size(), fp);
	}
	else
		gpi->desL.writeListBin(fp);

	imBN.grabaComprGris(fp, NULL, COMPR ? -1 : 0);

	dt=L_TIME()-t1; // Tiempo de calculo
	fwrite(&dt,sizeof(dt),1,fp);
	fwrite(&largoNomArch, sizeof(largoNomArch), 1, fp);
	fwrite(argv[1], 1, largoNomArch+1, fp);

	gMi.buscaCamposOrientacionFrecuencia();
	gMi.grabaCamposOrientacionFrecuencia(fp, COMPR);

	fclose(fp); // cierre Template

	fp=fopen(argv[4],"w"); // Output
	if (fp==NULL)
	{
		printf("No se puede generar archivo output en ENROLL\n");
		return 1;
	}
	t2=L_TIME();
	fprintf(fp, "%s %s XXX %.2fs OK\n", argv[1], argv[2], t2-t1);
	fclose(fp); // cierre Output

	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_enroll");
	return 0;
}

int L_ProgramasRecPuntInt::main_FVC_match(int argc, char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_match");

	if (argc!=5 && argc!=6)
	{
		printf(
			"\n"//0    1            2           3          4           5
			"  MATCH imagefile templatefile configfile outputfile [scorefile]\n"
			"    imagefile: una imagen de test\n"
			"    templatefile: archivo (datos de la imagen de referencia)\n"
			"    configfile: archivo.txt (parametros)\n"
			"    outputfile: archivo que indica la eficiencia del proceso, y si las imagenes calzan\n"
			"    scorefile: archivo opcional que guarda informacion de puntajes de la transformacion (si caractRes=true)\n"
			"Al usar usaMINEXT=true, se pueden agregar los archivos:\n"
			"    imagefile_min.txt    (minucias)\n"
			"    imagefile_min_OF.bmp (campos de orientacion y frecuencia en canales R, G)\n"
			"Esos archivos se usan para poder extraer caracteristicas extra\n"
			"Los parametros genImCalces, genImFlechasOriFrec, genListadoCalces, caractRes permiten generar archivos extra\n"
			"El parametro pausaAlFinal produce una pausa al final del proceso (presionar ENTER)\n"
			);
#ifndef __ATLIMAGE_H__
		printf("  --- Esta funcion read solo archivos BMP\n");
#else
		printf("  --- Esta version read varios formatos de imagen\n");
#endif
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_match");
		return 1;
	}

	L_ImageRGBUchar imRGB;
	L_ImageGrayUchar imBNRef, imCalidad;
	L_ImageGrayDouble imRef, imPru;

	L_FileName d(argv[1]);
	d.ext = "txt";
	d.name = d.name + "_min";

	// Variables usadas para el reconocedor
	L_GenDescrPuntInt *gpi;
	L_MinuciasExterno gMi(d.c_str());
#ifdef __L_LOWE_H_
	L_LoweGenDescr gLo;
#else
	L_PDoG_SIFT gLo;
#endif
	L_PHarris_SIFT gHa;
	L_Harris_SIFT gHa1;
	L_GenDescrArr gUnion(4);
	L_BusArr<L_KdTreeBBF> kdt(4); // 3 posibles puntos: DoG, minucias_termino, minucias_bifurcacion. Cada tipo build su propio KdTreeBBF
	L_GenTransfAfinHough hough;
	L_RecPuntInt rpi;
	L_ParamManagerLocal parLoc("FVC",5);
	L_ParamLabelList parL;
	L_DescriptorLista desL;
	L_FiltroDescriptoresGrupo filtGr;
	L_FiltroDesPolig filtPol(6);
	// Parametros para controlar los tipos de descriptores a usar, los calces y los dibujos generados
	bool v_activo=false;
	int v1_x=0, v1_y=0, v2_x=0, v2_y=0, v3_x=0, v3_y=0, v4_x=0, v4_y=0;
	bool genImCalces=false;
	bool genImFlechasOriFrec=false;
	bool genListadoCalces=false;
	bool caractRes=false;
	bool pausaAlFinal=true;
	bool usaSIFT=true;
	bool usaMINEXT=false;
	bool usaHARRIS=false;
	bool usaHARRIS1=false;
	// Variables usadas para comparar campos de orientacion y frecuencia
	double camposOF_minLongOnda=5;
	double camposOF_maxLongOnda=22;
	bool camposOF_v_activo=false;
	int camposOF_v1_x=0, camposOF_v1_y=0;
	int camposOF_v2_x=0, camposOF_v2_y=0;
	int camposOF_v3_x=0, camposOF_v3_y=0;
	int camposOF_v4_x=0, camposOF_v4_y=0;
	bool OFvalido;

	bool COMPR;
	FILE *fp;
	int largoNomArch;
	char nomImRef[1000];
	double t1, t2, dt;

	t1=L_TIME();
	gMi.agregaEspejo=true;

	// Configurar parametros
	parLoc.addFrom("v_activo",&v_activo);
	parLoc.addFrom("v1_x",&v1_x);
	parLoc.addFrom("v1_y",&v1_y);
	parLoc.addFrom("v2_x",&v2_x);
	parLoc.addFrom("v2_y",&v2_y);
	parLoc.addFrom("v3_x",&v3_x);
	parLoc.addFrom("v3_y",&v3_y);
	parLoc.addFrom("v4_x",&v4_x);
	parLoc.addFrom("v4_y",&v4_y);
	parLoc.addFrom("genImCalces",&genImCalces);
	parLoc.addFrom("genImFlechasOriFrec",&genImFlechasOriFrec);
	parLoc.addFrom("genListadoCalces",&genListadoCalces);
	parLoc.addFrom("caractRes",&caractRes);
	parLoc.addFrom("pausaAlFinal",&pausaAlFinal);
	parLoc.addFrom("usaMINEXT",&usaMINEXT);
	parLoc.addFrom("usaSIFT",&usaSIFT);
	parLoc.addFrom("usaHARRIS",&usaHARRIS);
	parLoc.addFrom("usaHARRIS1",&usaHARRIS1);
	// Parametros para comparacion de campos de orientacion y frecuencia
	parLoc.addFrom("camposOF_minLongOnda",&camposOF_minLongOnda);
	parLoc.addFrom("camposOF_maxLongOnda", &camposOF_maxLongOnda);
	parLoc.addFrom("camposOF_v_activo", &camposOF_v_activo);
	parLoc.addFrom("camposOF_v1_x", &camposOF_v1_x);
	parLoc.addFrom("camposOF_v1_y", &camposOF_v1_y);
	parLoc.addFrom("camposOF_v2_x", &camposOF_v2_x);
	parLoc.addFrom("camposOF_v2_y", &camposOF_v2_y);
	parLoc.addFrom("camposOF_v3_x", &camposOF_v3_x);
	parLoc.addFrom("camposOF_v3_y", &camposOF_v3_y);
	parLoc.addFrom("camposOF_v4_x", &camposOF_v4_x);
	parLoc.addFrom("camposOF_v4_y", &camposOF_v4_y);

	// Agregar sub-generadores de puntos a gUnion
	gUnion.agregaReferencia(&gMi);
	gUnion.agregaReferencia(&gLo);
	gUnion.agregaReferencia(&gHa);
	gUnion.agregaReferencia(&gHa1);
	gpi=&gUnion;
	rpi.fijaReconocedor(gpi, &kdt, &hough, &filtGr);
	parLoc.addChildren(&rpi.paramRPI);
	parL.readFile(fopen(argv[3],"r")); // configfile
	parLoc.updateValues_me_and_children(parL);
	// Leer todos los parametros y volverlos a write en el archivo
	parLoc.buildListOfValues_me_and_children(parL);
	parL.saveFile(fopen(argv[3],"w"));
	//activar sub-generadores de puntos
	gUnion.activo[0]=usaMINEXT;
	gUnion.activo[1]=usaSIFT;
	gUnion.activo[2]=usaHARRIS;
	gUnion.activo[3]=usaHARRIS1;

	if (v_activo) // Filtrar calces cuyos puntos salen de una cierta ventana (4 vertices)
	{
		filtGr.addFrom(&filtPol);
		filtPol.agregaVertice(v1_x,v1_y);
		filtPol.agregaVertice(v2_x,v2_y);
		filtPol.agregaVertice(v3_x,v3_y);
		filtPol.agregaVertice(v4_x,v4_y);
	}

	if (imRGB.readImage(argv[1])==0) // Imagefile
	{
		printf("No se puede abrir imagefile en MATCH\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_match");
		return 1;
	}

	imPru=imRGB; // imagefile

	fp=fopen(argv[2],"rb"); // templatefile
	if (fp==NULL)
	{
		printf("No se puede abrir archivo template en MATCH\n");
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_match");
		return 1;
	}
	// fp: desL, imBN, dt, largoNomArch, argv[1], camposOrientacionFrecuencia
	// Primero ver si esta comprimido
	COMPR = (fgetc(fp) != 0);
	if (COMPR)
	{

		std::vector<char> bufM, bufC;
		int enee;
		if (fread(&enee, sizeof(enee), 1, fp) != 1)
			{printf("main_FVC_match: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		bufC.resize(enee);
		if (fread(&(bufC[0]), sizeof(char), bufC.size(), fp) != bufC.size())
			{printf("main_FVC_match: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		L_HuffmanEncoder::decodeAll(bufC, bufM);
		desL.readListMem(bufM, 0, &imRef, gpi);
	}
	else
		desL.readListBin(fp, &imRef, gpi); // Se uso imPru como puntero-imagen para los descriptores del template. Cuidado con esta linea
	imBNRef.leeComprGris(fp);
	imRef = imBNRef;
	if (fread(&dt,sizeof(dt),1,fp) != 1)
		{printf("main_FVC_match: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	if (fread(&largoNomArch, sizeof(largoNomArch), 1, fp) != 1)
		{printf("main_FVC_match: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	if (fread(nomImRef, 1, L_MIN(largoNomArch+1, 999), fp) != largoNomArch+1)
		{printf("main_FVC_match: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	L_MinuciasExterno otherMi(NULL,NULL,NULL); // Solo para cargar los campos de frecuencia y orientacion
	otherMi.leeCamposOrientacionFrecuencia(fp);

	fclose(fp); // cierra templatefile

	OFvalido = gMi.buscaCamposOrientacionFrecuencia();

	nomImRef[999]=0; // Por si acaso
	imRef=imBNRef;
	rpi.agregaImagenEntrenCopiaObjetoDe(imRef, 0.5, &desL);

	rpi.finalizaEntren();
	rpi.procesaImagenesPrueba(imPru, 0.5);

	// Calculo de parecido de campos de orientacion y frecuencia
	double difOri=M_PI/2, difFrec=1.0, corrFrec=-1.0; // Valores muy malos
	int OF_vx[4], OF_vy[4];
	OF_vx[0]=camposOF_v1_x;	OF_vy[0]=camposOF_v1_y;
	OF_vx[1]=camposOF_v2_x;	OF_vy[1]=camposOF_v2_y;
	OF_vx[2]=camposOF_v3_x;	OF_vy[2]=camposOF_v3_y;
	OF_vx[3]=camposOF_v4_x;	OF_vy[3]=camposOF_v4_y;
	if (rpi.traL.size()>0)
		OFvalido = gMi.comparaCamposOrientacionFrec(difOri, difFrec, corrFrec, *(L_TransfAfinPI2D *)rpi.traL.root->c, otherMi, 1/camposOF_maxLongOnda, 1/camposOF_minLongOnda, camposOF_v_activo, OF_vx, OF_vy);
	else
		OFvalido = false;

	if (genImFlechasOriFrec && OFvalido)
	{
		gMi.grabaFlechasMinucias("flechasOriFrecRef.bmp", imRef);
		otherMi.grabaFlechasMinucias("flechasOriFrecPru.bmp", imPru);
	}

	if (genImCalces) // Esto es medio complicado cambiarlo... ademas no afecta al resto
	{
		L_ShapeArray lins;
		rpi.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
		L_Array<const L_ImageGrayDouble *> arr(2);
		arr[0] = &imPru;
		arr[1] = &imRef;
		imRGB.genDrawingMatches(arr, lins);
		imRGB.writeBMP("recon.bmp");
		lins.resize(0);
		rpi.calL.genLineas(lins);
		imRGB.genDrawingMatches(arr, lins);
		imRGB.writeBMP("todcal.bmp");
	}

	if (genListadoCalces) // Lista todos los calces de todas las transformaciones
	{
		FILE *archivo = fopen("listaCal.txt","w");
		if (archivo==NULL)
			printf("Error: no se pudo crear el archivo listaCal.txt\n");
		else
		{
			L_CalceNodo *cal;
			if (rpi.traL.size()>0)
			{
				fprintf(archivo, "%ld calces\n", rpi.traL.root->c->calL.size());
				for (cal = rpi.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
					fprintf(archivo, "%.2f %.2f %.2f %.2f\n", cal->c.dRef->x0, cal->c.dRef->y0, cal->c.dPru->x0, cal->c.dPru->y0);
			}
			fclose(archivo);
		}
	}

	fp=fopen(argv[4],"w"); // outputfile
	if (fp==NULL)
	{
		printf("No se pudo crear el archivo de output %s\n", argv[4]);
	}
	else
	{
		double punt = (rpi.traL.size()>0) ? rpi.traL.root->c->car.punt - rpi.traL.root->c->car.ransac : 0;
		t2=L_TIME();
		// imagefile templatefile matching_time result similarity
		fprintf(fp, "%s %s %.3f[s] %s %.4f\n", argv[1], argv[2], t2-t1, "OK", punt);
		fclose(fp);
	}

	if (caractRes || argc==6) // Escribir las caracteristicas extraidas
	{
		if (argc==6)
			fp=fopen(argv[5],"w");
		else
			fp=stdout;

		if (fp!=NULL)
		{
			// tiempo num_calces certeza n_trans votos prob ransac corrPix corrLin
			
			L_TransfAfinPI2D *tr;
			L_TransfAfinPI2D trVacio;
			L_String caractTxt;
			// Caracteristicas extraidas aca
			double tiempo;
			double n_transf;
			double n_cal;
			double escMax;
			double escMin;
			double rot;
			int i;

			tr = (rpi.traL.size()>0) ? (L_TransfAfinPI2D *)(rpi.traL.root->c) : &trVacio;
			caractTxt = tr->car.imprimeCaract();

			tiempo = dt + L_TIME()-t1;
			n_transf = rpi.traL.size();
			n_cal = rpi.caractRec.nCal;
			if (tr!=&trVacio)
				tr->descomponeSVD(escMax, escMin, rot);
			else
				escMax = escMin = rot = 0;
			// Datos

			fprintf(fp, "%.4f", n_cal); // Datos calculados aca

			for (i=0; i<tr->car.n; i++)
				fprintf(fp, " %.4f", tr->car.el(i)); // Datos de la transformacion principal

			fprintf(fp, " %.4f %.4f %.4f %.4f %.4f", escMax, escMin, rot, n_transf, tiempo); // Datos de la comparacion de imagenes de orientacion y frecuencia
			fprintf(fp, " %.4f %.4f %.4f", difOri, difFrec, corrFrec);
			fprintf(fp, "\n");
			// Texto descriptivo
			fprintf(fp, "n_cal ");
			fprintf(fp, "%s ", caractTxt.c_str());
			fprintf(fp, "escMax escMin rot n_transf tiempo");
			fprintf(fp, "difOri difFrec corrFrec");
			fprintf(fp, "\n");
			if (argc==6)
				fclose(fp);
		}
	}

	if (pausaAlFinal)
	{
		printf("Presione ENTER para terminar\n");
		getchar();
	}

	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_match");
	return 0;
}

int L_ProgramasRecPuntInt::main_FVC_matchT(int argc,char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");

	if (argc!=6)
	{
		printf(
			"\n" //0     1              2            3           4         5
			"  MATCHT templatefile list_template configfile outputfile scorefile\n\n"
			"    templatefile: un template de imagen de test\n"
			"    list_template: archivo de texto que indica los template de referencia a usar\n"
			"    configfile: archivo.txt (parametros)\n"
			"    outputfile: archivo que indica la eficiencia del proceso, y si las imagenes calzan\n"
			"scorefile: archivo opcional que guarda informacion de puntajes de la transformacion\n"
			"Los parametros genImCalces, genImFlechasOriFrec, genListadoCalces, caractRes permiten generar archivos extra\n"
			"El parametro pausaAlFinal produce una pausa al final del proceso (presionar ENTER)\n"
			);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");
		return 1;
	}

	L_Array<L_String> templRef;
	char linea[1000];
	bool COMPR;
	FILE *fp, *fpScoreFile;

	fp=fopen(argv[2], "r"); // list_template
	if (fp == NULL)
	{
		printf("No se encuentra el archivo list_template = %s\n", argv[2]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");
		return 2;
	}

	while (!feof(fp))
	{
		if (fgets(linea, 1000, fp) == NULL) // Fin de archivo
			break;
		templRef.resize(templRef.size()+1);
		templRef[templRef.size()-1] = linea;
		templRef[templRef.size()-1].erase_spacing_characters_at_beginning_and_at_ending();
		if (templRef[templRef.size()-1][0] == 0) // La linea contenia un name de archivo
			templRef.resize(templRef.size()-1);
	}
	fclose(fp); // cierra list_template

	// Aca templRef contiene los nombres de los templates que se deben read

	// Se van a usar los mismos generadores como referencia para todos los descriptores (referencia y test)
	L_GenDescrPuntInt *gpi;
	L_MinuciasExterno gMi(NULL, NULL, NULL); // No build los puntos, es solo para los parametros
#ifdef __L_LOWE_H_
	L_LoweGenDescr gLo;
#else
	L_PDoG_SIFT gLo;
#endif
	L_PHarris_SIFT gHa;
	L_Harris_SIFT gHa1;
	L_GenDescrArr gUnion(4);
	L_BusArr<L_KdTreeBBF> kdt(4); // 3 posibles puntos: DoG, minucias_termino, minucias_bifurcacion. Cada tipo build su propio KdTreeBBF
	L_GenTransfAfinHough hough;
	L_RecPuntInt rpi;
	L_ParamManagerLocal parLoc("FVC",5);
	L_ParamLabelList parL;
	L_FiltroDescriptoresGrupo filtGr;
	L_FiltroDesPolig filtPol(6);
	// Parametros para controlar los tipos de descriptores a usar, los calces y los dibujos generados
	bool v_activo=false;
	int v1_x=0, v1_y=0, v2_x=0, v2_y=0, v3_x=0, v3_y=0, v4_x=0, v4_y=0;
	bool genImCalces=false;
	bool genListadoCalces=false;
	bool usaSIFT=true;
	bool usaMINEXT=false;
	bool usaHARRIS=false;
	bool usaHARRIS1=false;
	// Variables usadas para comparar campos de orientacion y frecuencia
	double camposOF_minLongOnda=5;
	double camposOF_maxLongOnda=22;
	bool camposOF_v_activo=false;
	int camposOF_v1_x=0, camposOF_v1_y=0;
	int camposOF_v2_x=0, camposOF_v2_y=0;
	int camposOF_v3_x=0, camposOF_v3_y=0;
	int camposOF_v4_x=0, camposOF_v4_y=0;
	bool OFvalido;
	L_String caractTxt;

	int largoNomArch;
	char nomImRef[1000], nomImPru[1000];
	double t1, t2, dtPru, dtRef;

	t1=L_TIME();
	gMi.agregaEspejo=true;

	// Configurar parametros
	parLoc.addFrom("v_activo",&v_activo);
	parLoc.addFrom("v1_x",&v1_x);
	parLoc.addFrom("v1_y",&v1_y);
	parLoc.addFrom("v2_x",&v2_x);
	parLoc.addFrom("v2_y",&v2_y);
	parLoc.addFrom("v3_x",&v3_x);
	parLoc.addFrom("v3_y",&v3_y);
	parLoc.addFrom("v4_x",&v4_x);
	parLoc.addFrom("v4_y",&v4_y);
	parLoc.addFrom("genListadoCalces",&genListadoCalces);
	parLoc.addFrom("genImCalces",&genImCalces);
	parLoc.addFrom("usaMINEXT",&usaMINEXT);
	parLoc.addFrom("usaSIFT",&usaSIFT);
	parLoc.addFrom("usaHARRIS",&usaHARRIS);
	parLoc.addFrom("usaHARRIS1",&usaHARRIS1);
	// Parametros para comparacion de campos de orientacion y frecuencia
	parLoc.addFrom("camposOF_minLongOnda",&camposOF_minLongOnda);
	parLoc.addFrom("camposOF_maxLongOnda", &camposOF_maxLongOnda);
	parLoc.addFrom("camposOF_v_activo", &camposOF_v_activo);
	parLoc.addFrom("camposOF_v1_x", &camposOF_v1_x);
	parLoc.addFrom("camposOF_v1_y", &camposOF_v1_y);
	parLoc.addFrom("camposOF_v2_x", &camposOF_v2_x);
	parLoc.addFrom("camposOF_v2_y", &camposOF_v2_y);
	parLoc.addFrom("camposOF_v3_x", &camposOF_v3_x);
	parLoc.addFrom("camposOF_v3_y", &camposOF_v3_y);
	parLoc.addFrom("camposOF_v4_x", &camposOF_v4_x);
	parLoc.addFrom("camposOF_v4_y", &camposOF_v4_y);

	// Agregar sub-generadores de puntos a gUnion
	gUnion.agregaReferencia(&gMi);
	gUnion.agregaReferencia(&gLo);
	gUnion.agregaReferencia(&gHa);
	gUnion.agregaReferencia(&gHa1);
	gpi=&gUnion;
	rpi.fijaReconocedor(gpi, &kdt, &hough, &filtGr);
	parLoc.addChildren(&rpi.paramRPI);
	parL.readFile(fopen(argv[3],"r")); // configfile
	parLoc.updateValues_me_and_children(parL);
	// Leer todos los parametros y volverlos a write en el archivo
	parLoc.buildListOfValues_me_and_children(parL);
	parL.saveFile(fopen(argv[3],"w"));
	//activar sub-generadores de puntos
	gUnion.activo[0]=usaMINEXT;
	gUnion.activo[1]=usaSIFT;
	gUnion.activo[2]=usaHARRIS;
	gUnion.activo[3]=usaHARRIS1;

	if (v_activo) // Filtrar calces cuyos puntos salen de una cierta ventana (4 vertices)
	{
		filtGr.addFrom(&filtPol);
		filtPol.agregaVertice(v1_x,v1_y);
		filtPol.agregaVertice(v2_x,v2_y);
		filtPol.agregaVertice(v3_x,v3_y);
		filtPol.agregaVertice(v4_x,v4_y);
	}

	// Leer template de test
	L_DescriptorLista desLPru, desLRef, desLPruCopia;
	L_ImageGrayUchar imBNPru, imBNRef;
	L_ImageGrayDouble imPru, imRef;

	fp=fopen(argv[1],"rb"); // templatefile
	if (fp==NULL)
	{
		printf("No se puede abrir archivo template %s en MATCHT\n", argv[1]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");
		return 1;
	}
	// fp: desL, imBN, dt, largoNomArch, argv[1], camposOrientacionFrecuencia
	COMPR = (fgetc(fp) != 0);
	if (COMPR)
	{

		std::vector<char> bufM, bufC;
		int enee;
		if (fread(&enee, sizeof(enee), 1, fp) != 1)
			{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		bufC.resize(enee);
		if (fread((&bufC[0]), sizeof(char), bufC.size(), fp) != bufC.size())
			{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		L_HuffmanEncoder::decodeAll(bufC, bufM);
		desLPru.readListMem(bufM, 0, &imPru, gpi);
	}
	else
		desLPru.readListBin(fp, &imPru, gpi); // Se uso imPru como puntero-imagen para los descriptores del template. Cuidado con esta linea

	imBNPru.leeComprGris(fp);
	imPru = imBNPru;
	if (fread(&dtPru,sizeof(dtPru),1,fp) != 1)
		{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	if (fread(&largoNomArch, sizeof(largoNomArch), 1, fp) != 1)
		{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	if (fread(nomImPru, 1, L_MIN(largoNomArch+1, 999), fp) != largoNomArch+1)
		{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
	L_MinuciasExterno otherMi(NULL,NULL,NULL); // Solo para cargar los campos de frecuencia y orientacion
	gMi.leeCamposOrientacionFrecuencia(fp);
	fclose(fp); // cierra templatefile

	fpScoreFile = fopen(argv[5], "w");
	if (fpScoreFile==NULL)
	{
		printf("No se puede abrir archivo score %s en MATCHT\n", argv[5]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");
		return 1;
	}

	// Parametros para parecido de campos de orientacion y frecuencia
	double difOri=M_PI/2, difFrec=1.0, corrFrec=-1.0; // Valores muy malos
	int OF_vx[4], OF_vy[4];
	OF_vx[0]=camposOF_v1_x;	OF_vy[0]=camposOF_v1_y;
	OF_vx[1]=camposOF_v2_x;	OF_vy[1]=camposOF_v2_y;
	OF_vx[2]=camposOF_v3_x;	OF_vy[2]=camposOF_v3_y;
	OF_vx[3]=camposOF_v4_x;	OF_vy[3]=camposOF_v4_y;

	int iTemplate;

	// Ahora a read cada template de referencia en list_template
	for (iTemplate = 0; iTemplate < (int)templRef.size(); iTemplate++)
	{
		double t1 = L_TIME();
		fp = fopen(templRef[iTemplate].c_str(), "rb");
		if (fp==NULL)
		{
			printf("No se puede abrir archivo template %s en MATCHT\n", templRef[iTemplate].c_str());
			L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");
			return 1;
		}
		COMPR = (fgetc(fp) != 0);
		if (COMPR)
		{

			std::vector<char> bufM, bufC;
			int enee;
			if (fread(&enee, sizeof(enee), 1, fp) != 1)
				{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
			bufC.resize(enee);
			if (fread(&(bufC[0]), sizeof(char), bufC.size(), fp) != bufC.size())
				{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
			L_HuffmanEncoder::decodeAll(bufC, bufM);
			desLRef.readListMem(bufM, 0, &imRef, gpi);
		}
		else
			desLRef.readListBin(fp, &imRef, gpi); // Se uso imPru como puntero-imagen para los descriptores del template. Cuidado con esta linea
		desLRef.fijaNumObj(0);
		imBNRef.leeComprGris(fp);
		imRef = imBNRef;
		if (fread(&dtRef,sizeof(dtRef),1,fp) != 1)
			{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		if (fread(&largoNomArch, sizeof(largoNomArch), 1, fp) != 1)
			{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		if (fread(nomImRef, 1, L_MIN(largoNomArch+1, 999), fp) != largoNomArch+1)
			{printf("main_FVC_matchT: Formato de archivo incorrecto\n"); fclose(fp); return 1;}
		L_MinuciasExterno otherMi(NULL,NULL,NULL); // Solo para cargar los campos de frecuencia y orientacion
		otherMi.leeCamposOrientacionFrecuencia(fp);
		fclose(fp); // cierra templatefile

		// Ahora el reconocimiento

		desLPruCopia.clear();
		desLPru.copyListOn(desLPruCopia);

		rpi.agregaImagenesEntren_usa_desL(desLRef, &imRef);
		rpi.finalizaEntren();
		rpi.procesaImagenesPrueba_usa_desL(desLPruCopia, &imPru);

		if (rpi.traL.size()>0)
			OFvalido = gMi.comparaCamposOrientacionFrec(difOri, difFrec, corrFrec, *(L_TransfAfinPI2D *)rpi.traL.root->c, otherMi, 1/camposOF_maxLongOnda, 1/camposOF_minLongOnda, camposOF_v_activo, OF_vx, OF_vy);
		else
			OFvalido = false;

		if (genImCalces) // Esto es medio complicado cambiarlo... ademas no afecta al resto
		{
			L_ImageRGBUchar imRGB;
			L_FileName ar1, ar2;
			ar1.asign(templRef[iTemplate].c_str());
			ar1.name = ar1.name + "_tra";
			ar1.ext = "bmp";
			ar2.asign(templRef[iTemplate].c_str());
			ar2.name = ar2.name + "_cal";
			ar2.ext = "bmp";
			L_ShapeArray lins;
			rpi.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
			L_Array<const L_ImageGrayDouble *> arr(2);
			arr[0] =  &imPru;
			arr[1] =  &imRef;
			imRGB.genDrawingMatches(arr, lins);
			imRGB.writeBMP(ar1.dirnameext());
			lins.resize(0);
			rpi.calL.genLineas(lins);
			imRGB.genDrawingMatches(arr, lins);
			imRGB.writeBMP(ar2.dirnameext());
		}

		if (genListadoCalces) // Lista todos los calces de todas las transformaciones
		{
			L_FileName ar;
			ar.asign(templRef[iTemplate].c_str());
			ar.name = ar.name + "_cal";
			ar.ext = "txt";
			FILE *archivo = fopen(ar.dirnameext(),"w");
			if (archivo==NULL)
				printf("Error: no se pudo crear el archivo %s\n", ar.dirnameext());
			else
			{
				L_CalceNodo *cal;
				if (rpi.traL.size()>0)
				{
					fprintf(archivo, "%ld calces\n", long(rpi.traL.root->c->calL.size()));
					for (cal = rpi.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
						fprintf(archivo, "%.2f %.2f %.2f %.2f\n", cal->c.dRef->x0, cal->c.dRef->y0, cal->c.dPru->x0, cal->c.dPru->y0);
				}
				fclose(archivo);
			}
		}

		fp=fopen(argv[4],"w"); // outputfile
		if (fp==NULL)
		{
			printf("No se pudo crear el archivo de output %s\n", argv[4]);
		}
		else
		{
			double punt = (rpi.traL.size()>0) ? rpi.traL.root->c->car.punt - rpi.traL.root->c->car.ransac : 0;
			t2=L_TIME();
			// imagefile templatefile matching_time result similarity
			fprintf(fp, "%s %s %.3f[s] %s %.4f\n", argv[1], argv[2], t2-t1, "OK", punt);
			fclose(fp);
		}

		// Escribir las caracteristicas extraidas en la deteccion
		// tiempo num_calces certeza n_trans votos prob ransac corrPix corrLin
		L_TransfAfinPI2D *tr;
		L_TransfAfinPI2D trVacio;
		// Caracteristicas extraidas aca
		double tiempo;
		double n_transf;
		double n_cal;
		double escMax;
		double escMin;
		double rot;
		int i;

		tr = (rpi.traL.size()>0) ? (L_TransfAfinPI2D *)(rpi.traL.root->c) : &trVacio;

		n_transf = rpi.traL.size();
		n_cal = rpi.caractRec.nCal;
		if (tr!=&trVacio)
			tr->descomponeSVD(escMax, escMin, rot);
		else
			escMax = escMin = rot = 0;
		tiempo = dtPru + dtRef + L_TIME()-t1;  // Se le suman los tiempos guardados en los templates


		fprintf(fpScoreFile, "%.4f", n_cal); // Datos calculados aca

		for (i=0; i<tr->car.n; i++)
			fprintf(fpScoreFile, " %.4f", tr->car.el(i)); // Datos de la transformacion principal

		fprintf(fpScoreFile, " %.4f %.4f %.4f %.4f %.4f", escMax, escMin, rot, n_transf, tiempo); // Datos de la comparacion de imagenes de orientacion y frecuencia
		fprintf(fpScoreFile, " %.4f %.4f %.4f", difOri, difFrec, corrFrec);
		fprintf(fpScoreFile, "\n");

		rpi.reseteaPrueba();
		rpi.reseteaEntren();
	}
	// Texto descriptivo
	caractTxt = L_Transf2D_caract::imprimeCaract();
	fprintf(fpScoreFile, "n_cal ");
	fprintf(fpScoreFile, "%s", caractTxt.c_str());
	fprintf(fpScoreFile, "escMax escMin rot n_transf tiempo ");
	fprintf(fpScoreFile, "difOri difFrec corrFrec");
	fprintf(fpScoreFile, "\n");
	fclose(fpScoreFile);

	//rpi.ti.escribeResumen(stdout);

	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT");

	return 0;
}



int L_ProgramasRecPuntInt::main_FVC_matchT2(int argc,char *argv[])
{
	L_PUSH_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");

	if (argc!=6)
	{
		printf(
			"\n" //0         1                 2            3           4         5
			"  MATCHT2 list_template_test list_template configfile outputfile scorefile\n\n"
			"    list_template_test: archivo de texto que indica los template de test\n"
			"    list_template: archivo de texto que indica los template de referencia a usar\n"
			"    configfile: archivo.txt (parametros)\n"
			"    outputfile: archivo que indica la eficiencia del proceso, y si las imagenes calzan\n"
			"scorefile: archivo opcional que guarda informacion de puntajes de la transformacion\n"
			"Los parametros genImCalces, genImFlechasOriFrec, genListadoCalces, caractRes permiten generar archivos extra\n"
			"El parametro pausaAlFinal produce una pausa al final del proceso (presionar ENTER)\n"
			);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
		return 1;
	}

	L_Array<L_String> templPru, templRef;
	char linea[1000];
	FILE *fp, *fpScoreFile;

	////////
	// Lectura lista de templates de test
	fp=fopen(argv[1], "r"); // list_template
	if (fp == NULL)
	{
		printf("No se encuentra el archivo list_template = %s\n", argv[2]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
		return 1;
	}

	while (!feof(fp))
	{
		if (fgets(linea, 1000, fp) == NULL) // Fin de archivo
			break;
		templPru.resize(templPru.size()+1);
		templPru[templPru.size()-1] = linea;
		templPru[templPru.size()-1].erase_spacing_characters_at_beginning_and_at_ending();
		if (templPru[templPru.size()-1][0] == 0) // La linea contenia un name de archivo
			templPru.resize(templPru.size()-1);
	}
	fclose(fp); // cierra list_template

	//////
	// Lectura lista de templates de referencia
	fp=fopen(argv[2], "r"); // list_template
	if (fp == NULL)
	{
		printf("No se encuentra el archivo list_template = %s\n", argv[2]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
		return 2;
	}

	while (!feof(fp))
	{
		if (fgets(linea, 1000, fp) == NULL) // Fin de archivo
			break;
		templRef.resize(templRef.size()+1);
		templRef[templRef.size()-1] = linea;
		templRef[templRef.size()-1].erase_spacing_characters_at_beginning_and_at_ending();
		if (templRef[templRef.size()-1][0] == 0) // La linea contenia un name de archivo
			templRef.resize(templRef.size()-1);
	}
	fclose(fp); // cierra list_template

	// Aca (templPru, templRef) contienen los nombres de los templates que se deben read

	//////
	//  Objetos a usar
	//  Se van a usar los mismos generadores como referencia para todos los descriptores (referencia y test)
	L_TemplateDescr tPru, tRef;
	L_BusArr<L_KdTreeBBF> kdt(4); // 3 posibles puntos: DoG, minucias_termino, minucias_bifurcacion. Cada tipo build su propio KdTreeBBF
	L_GenTransfAfinHough hough;
	L_RecPuntInt rpi;
	L_ParamManagerLocal parLoc("FVC",5);
	L_ParamLabelList parL;
	L_FiltroDescriptoresGrupo filtGr;
	L_FiltroDesPolig filtPol(6);
	// Parametros para controlar los tipos de descriptores a usar, los calces y los dibujos generados
	bool usaSIFT=true;
	bool usaMINEXT=false;
	bool usaHARRIS=false;
	bool usaHARRIS1=false;
	bool genImCalces=false;
	bool genListadoCalces=false;
	bool v_activo=false;
	int v1_x=0, v1_y=0, v2_x=0, v2_y=0, v3_x=0, v3_y=0, v4_x=0, v4_y=0;
	// Variables usadas para comparar campos de orientacion y frecuencia
	double camposOF_minLongOnda=5;
	double camposOF_maxLongOnda=22;
	bool camposOF_v_activo=false;
	int camposOF_v1_x=0, camposOF_v1_y=0;
	int camposOF_v2_x=0, camposOF_v2_y=0;
	int camposOF_v3_x=0, camposOF_v3_y=0;
	int camposOF_v4_x=0, camposOF_v4_y=0;
	bool OFvalido;
	L_String caractTxt;

	double t1, t2;

	// Configurar parametros
	parLoc.addFrom("v_activo",&v_activo);
	parLoc.addFrom("v1_x",&v1_x);
	parLoc.addFrom("v1_y",&v1_y);
	parLoc.addFrom("v2_x",&v2_x);
	parLoc.addFrom("v2_y",&v2_y);
	parLoc.addFrom("v3_x",&v3_x);
	parLoc.addFrom("v3_y",&v3_y);
	parLoc.addFrom("v4_x",&v4_x);
	parLoc.addFrom("v4_y",&v4_y);
	parLoc.addFrom("genListadoCalces",&genListadoCalces);
	parLoc.addFrom("genImCalces",&genImCalces);
	parLoc.addFrom("usaMINEXT",&usaMINEXT);
	parLoc.addFrom("usaSIFT",&usaSIFT);
	parLoc.addFrom("usaHARRIS",&usaHARRIS);
	parLoc.addFrom("usaHARRIS1",&usaHARRIS1);
	// Parametros para comparacion de campos de orientacion y frecuencia
	parLoc.addFrom("camposOF_minLongOnda",&camposOF_minLongOnda);
	parLoc.addFrom("camposOF_maxLongOnda", &camposOF_maxLongOnda);
	parLoc.addFrom("camposOF_v_activo", &camposOF_v_activo);
	parLoc.addFrom("camposOF_v1_x", &camposOF_v1_x);
	parLoc.addFrom("camposOF_v1_y", &camposOF_v1_y);
	parLoc.addFrom("camposOF_v2_x", &camposOF_v2_x);
	parLoc.addFrom("camposOF_v2_y", &camposOF_v2_y);
	parLoc.addFrom("camposOF_v3_x", &camposOF_v3_x);
	parLoc.addFrom("camposOF_v3_y", &camposOF_v3_y);
	parLoc.addFrom("camposOF_v4_x", &camposOF_v4_x);
	parLoc.addFrom("camposOF_v4_y", &camposOF_v4_y);

	// Agregar sub-generadores de puntos a gUnion
	tRef.fijaGeneradores(usaMINEXT, usaSIFT, usaHARRIS, usaHARRIS1);
	tPru.fijaGeneradores(usaMINEXT, usaSIFT, usaHARRIS, usaHARRIS1);
	rpi.fijaReconocedor(&tRef, &kdt, &hough, &filtGr);
	parLoc.addChildren(&rpi.paramRPI);
	parL.readFile(fopen(argv[3],"r")); // configfile
	L_ParamLabelList parL2;
	parL.copyListOn(parL2);
	parLoc.updateValues_me_and_children(parL);
	// Actualizar valores del template
	tPru.paramTempl.updateValues_me_and_children(parL2);
	tPru.pideParams()->updateValues_me_and_children(parL2);

	// Leer todos los parametros y volverlos a write en el archivo
	parLoc.buildListOfValues_me_and_children(parL);
	parL.saveFile(fopen(argv[3],"w"));

	if (v_activo) // Filtrar calces cuyos puntos salen de una cierta ventana (4 vertices)
	{
		filtGr.addFrom(&filtPol);
		filtPol.agregaVertice(v1_x,v1_y);
		filtPol.agregaVertice(v2_x,v2_y);
		filtPol.agregaVertice(v3_x,v3_y);
		filtPol.agregaVertice(v4_x,v4_y);
	}

	// Parametros para parecido de campos de orientacion y frecuencia
	double difOri=M_PI/2, difFrec=1.0, corrFrec=-1.0; // Valores muy malos
	int OF_vx[4], OF_vy[4];
	OF_vx[0]=camposOF_v1_x;	OF_vy[0]=camposOF_v1_y;
	OF_vx[1]=camposOF_v2_x;	OF_vy[1]=camposOF_v2_y;
	OF_vx[2]=camposOF_v3_x;	OF_vy[2]=camposOF_v3_y;
	OF_vx[3]=camposOF_v4_x;	OF_vy[3]=camposOF_v4_y;

	fpScoreFile = fopen(argv[5], "w");
	if (fpScoreFile==NULL)
	{
		printf("No se puede abrir archivo score %s en MATCHT\n", argv[5]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
		return 3;
	}

	// Texto descriptivo
	fprintf(fpScoreFile, "Archivo de test: %s\n", argv[1]);
	fprintf(fpScoreFile, "Archivo de referencia: %s\n", argv[2]);
	fprintf(fpScoreFile, "nImPru nImRef ");
	caractTxt = L_Transf2D_caract::imprimeCaract();
	fprintf(fpScoreFile, "n_cal ");
	fprintf(fpScoreFile, "%s", caractTxt.c_str());
	fprintf(fpScoreFile, "escMax escMin rot n_transf tiempo ");
	fprintf(fpScoreFile, "difOri difFrec corrFrec");
	fprintf(fpScoreFile, "\n");

	fp=fopen(argv[4],"w"); // outputfile
	if (fp==NULL)
	{
		printf("No se pudo crear el archivo de output %s\n", argv[4]);
		L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
		return 4;
	}

	tPru.nobj = -1;
	tRef.nobj = 0;

	// Ahora a read cada template de referencia en list_template
	int iTemplate, iTest;
	for (iTemplate = 0; iTemplate < (int)templRef.size(); iTemplate++)
	{
		if (tRef.leeTemplate(templRef[iTemplate].c_str()) == false)
		{
			L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
			return 5;
		}
		// Ahora el reconocimiento
		rpi.agregaImagenesEntren_usa_desL(tRef.desL, &tRef.im);
		rpi.finalizaEntren();
		for (iTest = 0; iTest < (int)templPru.size(); iTest++)
		{
			// Leer template de test
			t1=L_TIME();
			if (tPru.leeTemplate(templPru[iTest].c_str()) == false)
			{

				L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");
				return 6;
			}

			rpi.procesaImagenesPrueba_usa_desL(tPru.desL, &tPru.im);

			if (rpi.traL.size()>0)
				OFvalido = tPru.gMi.comparaCamposOrientacionFrec(difOri, difFrec, corrFrec, *(L_TransfAfinPI2D *)rpi.traL.root->c, tRef.gMi, 1/camposOF_maxLongOnda, 1/camposOF_minLongOnda, camposOF_v_activo, OF_vx, OF_vy);
			else
				OFvalido = false;

			if (genImCalces)
			{
				L_ImageRGBUchar imRGB;
				L_FileName ar1, ar2;
				ar1.asign(templRef[iTemplate].c_str());
				ar1.name = ar1.name + "_tra";
				ar1.ext = "bmp";
				ar2.asign(templRef[iTemplate].c_str());
				ar2.name = ar2.name + "_cal";
				ar2.ext = "bmp";
				rpi.dibujaTransformaciones(ar1.c_str(), 1, 1, 1, 1);
				rpi.dibujaCalcesIniciales(ar2.c_str());
			}

			if (genListadoCalces) // Lista todos los calces de todas las transformaciones
			{
				L_FileName ar;
				ar.asign(templRef[iTemplate].c_str());
				ar.name = ar.name + "_cal";
				ar.ext = "txt";
				FILE *archivo = fopen(ar.dirnameext(),"w");
				if (archivo==NULL)
					printf("Error: no se pudo crear el archivo %s\n", ar.dirnameext());
				else
				{
					L_CalceNodo *cal;
					if (rpi.traL.size()>0)
					{
						fprintf(archivo, "%ld calces\n", rpi.traL.root->c->calL.size());
						for (cal = rpi.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
							fprintf(archivo, "%.2f %.2f %.2f %.2f\n", cal->c.dRef->x0, cal->c.dRef->y0, cal->c.dPru->x0, cal->c.dPru->y0);
					}
					fclose(archivo);
				}
			}

			if (fp != NULL)
			{
				double punt = (rpi.traL.size()>0) ? rpi.traL.root->c->car.punt - rpi.traL.root->c->car.ransac : 0;
				t2=L_TIME();
				// imagefile templatefile matching_time result similarity
				fprintf(fp, "%s %s %.3f[s] %s %.4f\n", argv[1], argv[2], t2-t1, "OK", punt);
			}

			// Escribir las caracteristicas extraidas en la deteccion
			// tiempo num_calces certeza n_trans votos prob ransac corrPix corrLin
			L_TransfAfinPI2D *tr;
			L_TransfAfinPI2D trVacio;
			// Caracteristicas extraidas aca
			double tiempo;
			double n_transf;
			double n_cal;
			double escMax;
			double escMin;
			double rot;
			int i;

			tr = (rpi.traL.size()>0) ? (L_TransfAfinPI2D *)(rpi.traL.root->c) : &trVacio;

			n_transf = rpi.traL.size();
			n_cal = rpi.caractRec.nCal;
			if (tr!=&trVacio)
				tr->descomponeSVD(escMax, escMin, rot);
			else
				escMax = escMin = rot = 0;
			tiempo = tPru.dt + tRef.dt + L_TIME()-t1;  // Se le suman los tiempos guardados en los templates

			fprintf(fpScoreFile, "%d %d ", iTest+1, iTemplate+1);
			fprintf(fpScoreFile, "%.4f", n_cal); // Datos calculados aca

			for (i=0; i<tr->car.n; i++)
				fprintf(fpScoreFile, " %.4f", tr->car.el(i)); // Datos de la transformacion principal

			fprintf(fpScoreFile, " %.4f %.4f %.4f %.4f %.4f", escMax, escMin, rot, n_transf, tiempo); // Datos de la comparacion de imagenes de orientacion y frecuencia
			fprintf(fpScoreFile, " %.4f %.4f %.4f", difOri, difFrec, corrFrec);
			fprintf(fpScoreFile, "\n");

			rpi.reseteaPrueba();
		}
		rpi.reseteaEntren();
	}
	fclose(fpScoreFile);
	fclose(fp);

	//rpi.ti.escribeResumen(stdout);
	printf("Mean enroll time: %f\n", rpi.ti.gpiTime.getMean() + rpi.ti.gbdTime.getMean());
	printf("Mean match time: %f\n", rpi.ti.calTime.getMean() + rpi.ti.gtcTime.getMean());

	L_POP_EXECUTING_FN("L_ProgramasRecPuntInt::main_FVC_matchT2");

	return 0;
}

int L_ProgramasRecPuntInt::main_pruebaCruzada(int argc, char *argv[])
{
	if (argc!=5)
	{
		printf(
			"\n"
			"  PRUCRUZ imEntren.txt imPrueba.txt param.txt output.txt\n"
			"\n"
			"imEntren: archivo de texto con los nombres de las imagenes de entrenamiento\n"
			"imPrueba: archivo de texto con los nombres de las imagenes de test\n"
			"param.txt: archivo con los parametros del sistema\n"
			"output.txt: archivo con la output del sistema\n");
		return 0;
	}
	FILE *fpEntren, *fpPrueba, *fpSal;
	L_ImageGrayDouble *imEntren, *imPrueba;
	L_DescriptorLista *desL, lista;
	int nEntren, nPrueba, len, i, j;
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	L_KdTreeBBF tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	char str[500];
	double prob;
	int nVotos;
	int nTransf;

	// Contar cuantas imagenes de entrenamiento hay, luego leerlas
	fpEntren=fopen(argv[1],"r");
	if (fpEntren==NULL)
	{
		printf("No se encuentra el archivo \"%s\"\n", argv[1]);
		return 1;
	}
	nEntren=0;
	while(true)
	{
		if (fgets(str, 499, fpEntren) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		nEntren++;
	}
	rewind(fpEntren);
	imEntren=new L_ImageGrayDouble[nEntren];
	nEntren=0;
	while(true)
	{
		if (fgets(str, 499, fpEntren) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		len=(int)strlen(str);
		while (str[len-1]=='\n' || str[len-1]==' ')
			str[len-- - 1] = 0;
		if ( imEntren[nEntren].readImage(str) == 0 )
		{
			printf("No se encuentra la imagen de entrenamiento %s\n", str);
			return 2;
		}
		nEntren++;
	}
	fclose(fpEntren);

	// Contar cuantas imagenes de test hay, luego leerlas
	fpPrueba=fopen(argv[2],"r");
	if (fpPrueba==NULL)
	{
		printf("No se encuentra el archivo \"%s\"\n", argv[2]);
		return 1;
	}
	nPrueba=0;
	while(true)
	{
		if (fgets(str, 499, fpPrueba) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		nPrueba++;
	}
	rewind(fpPrueba);
	imPrueba=new L_ImageGrayDouble[nPrueba];
	nPrueba=0;
	while(true)
	{
		if (fgets(str, 499, fpPrueba) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		len=(int)strlen(str);
		while (str[len-1]=='\n' || str[len-1]==' ')
			str[len-- - 1] = 0;
		if ( imPrueba[nPrueba].readImage(str) == 0 )
		{
			printf("No se encuentra la imagen de test %s\n", str);
			return 2;
		}
		nPrueba++;
	}
	fclose(fpPrueba);

	// Formar el reconocedor
	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.leeParametrosPropaga(argv[3]);
	rec.grabaParametrosPropaga(argv[3]);
	desL = new L_DescriptorLista[nPrueba];
	sift.nobj=-1;

	// Calcular descriptores de test
	for (i=0; i<nPrueba; i++)
	{
		printf("Procesando pru %d (de %d)\n", i, nPrueba);
		sift.calcDescrPuntInt(imPrueba[i]);
		sift.desL.moveListTo(desL[i]);
	}

	// Abrir archivo de output
	fpSal = fopen(argv[4],"w");
	if (fpSal==NULL)
	{
		printf("No se puede abrir el archivo de output %s\n", argv[4]);
		return 3;
	}

	// ciclo de entrenamiento
	fprintf(fpSal, "Formato:  nEntren(0-...)  nPrueba(0-...)  prob  nVotos  nTransf\n");
	for (i=0; i<nEntren; i++)
	{
		rec.agregaImagenEntrenRobaObjetoDe(imEntren[i]);
		rec.finalizaEntren();
		// ciclo de reconocimiento
		for (j=0; j<nPrueba; j++)
		{
			printf("Entren: %d(%d) test: %d(%d)\n", i, nEntren, j, nPrueba);
			//desL[j].cloneListOn(lista);
			desL[j].copyListOn(lista);
			rec.procesaImagenesPrueba_usa_desL(lista, &imPrueba[j]);
			// Grabar resultados del reconocimiento
			if (rec.traL.root!=NULL)
			{
				prob=rec.traL.root->c->car.punt;
				nVotos=(int)rec.traL.root->c->calL.size();
			}
			else
			{
				prob=0;
				nVotos=0;
			}
			nTransf=(int)rec.traL.size();

			fprintf(fpSal, "%d %d %f %d %d\n", i, j, prob, nVotos, nTransf);
			rec.reseteaPrueba();
		}
		rec.reseteaEntren();
	}
	fclose(fpSal);
	// Destructores
	delete[] imEntren;
	delete[] imPrueba;
	delete[] desL;
	return 0;
}

int L_ProgramasRecPuntInt::main_pruebaCruzada_1v1(int argc, char *argv[])
{
	if (argc!=5)
	{
		printf(
			"\n"
			"  PRUCRUZ imEntren.txt imPrueba.txt param.txt output.txt\n"
			"\n"
			"imEntren: archivo de texto con los nombres de las imagenes de entrenamiento\n"
			"imPrueba: archivo de texto con los nombres de las imagenes de test\n"
			"param.txt: archivo con los parametros del sistema\n"
			"output.txt: archivo con la output del sistema\n");
		return 0;
	}
	FILE *fpEntren, *fpPrueba, *fpSal;
	L_ImageGrayDouble *imEntren, *imPrueba;
	int nEntren, nPrueba, len, i, j;
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	L_KdTreeBBF tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	char str[500];
	double prob;
	int nVotos;
	int nTransf;
	double t1, t2;

	// Contar cuantas imagenes de entrenamiento hay, luego leerlas
	fpEntren=fopen(argv[1],"r");
	if (fpEntren==NULL)
	{
		printf("No se encuentra el archivo \"%s\"\n", argv[1]);
		return 1;
	}
	nEntren=0;
	while(true)
	{
		if (fgets(str, 499, fpEntren) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		nEntren++;
	}
	rewind(fpEntren);
	imEntren=new L_ImageGrayDouble[nEntren];
	nEntren=0;
	while(true)
	{
		if (fgets(str, 499, fpEntren) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		len=(int)strlen(str);
		while (str[len-1]=='\n' || str[len-1]==' ')
			str[len-- - 1] = 0;
		if ( imEntren[nEntren].readImage(str) == 0 )
		{
			printf("No se encuentra la imagen de entrenamiento %s\n", str);
			return 2;
		}
		nEntren++;
	}
	fclose(fpEntren);

	// Contar cuantas imagenes de test hay, luego leerlas
	fpPrueba=fopen(argv[2],"r");
	if (fpPrueba==NULL)
	{
		printf("No se encuentra el archivo \"%s\"\n", argv[2]);
		return 1;
	}
	nPrueba=0;
	while(true)
	{
		if (fgets(str, 499, fpPrueba) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		nPrueba++;
	}
	rewind(fpPrueba);
	imPrueba=new L_ImageGrayDouble[nPrueba];
	nPrueba=0;
	while(true)
	{
		if (fgets(str, 499, fpPrueba) == NULL)
			break;
		if (*str==' ' || *str=='\n')
			continue;
		len=(int)strlen(str);
		while (str[len-1]=='\n' || str[len-1]==' ')
			str[len-- - 1] = 0;
		if ( imPrueba[nPrueba].readImage(str) == 0 )
		{
			printf("No se encuentra la imagen de test %s\n", str);
			return 2;
		}
		nPrueba++;
	}
	fclose(fpPrueba);

	// Comenzar a contar el tiempo
	t1 = L_TIME();

	// Formar el reconocedor
	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.leeParametrosPropaga(argv[3]);
	rec.grabaParametrosPropaga(argv[3]);

	// Abrir archivo de output
	fpSal = fopen(argv[4],"w");
	if (fpSal==NULL)
	{
		printf("No se puede abrir el archivo de output %s\n", argv[4]);
		return 3;
	}

	// ciclo de entrenamiento
	fprintf(fpSal, "Formato:  nEntren(0-...)  nPrueba(0-...)  prob  nVotos  nTransf\n");
	for (i=0; i<nEntren; i++)
	{
		rec.agregaImagenEntrenRobaObjetoDe(imEntren[i]);
		rec.finalizaEntren();
		// ciclo de reconocimiento
		for (j=0; j<nPrueba; j++)
		{
			printf("Entren: %d(%d) test: %d(%d)\n", i, nEntren, j, nPrueba);
			rec.procesaImagenesPrueba(imPrueba[j]);
			// Grabar resultados del reconocimiento
			if (rec.traL.root!=NULL)
			{
				prob=rec.traL.root->c->car.punt;
				nVotos=(int)rec.traL.root->c->calL.size();
			}
			else
			{
				prob=0;
				nVotos=0;
			}
			nTransf=(int)rec.traL.size();

			fprintf(fpSal, "%d %d %f %d %d\n", i, j, prob, nVotos, nTransf);
			rec.reseteaPrueba();
		}
		rec.reseteaEntren();
	}
	// terminar de contar el tiempo
	t2 = L_TIME();

	// mostrar el tiempo en pantalla e imprimir al archivo
	printf("Tiempo total: %.6f[s]\n", t2-t1);
	printf("Tiempo por cada par: %.6f[s]\n", (t2-t1)/(nEntren*(double)nPrueba));
	printf("Tiempo optimizado por cada par: %.6f[s]\n", rec.ti.calTime.getMean()+rec.ti.gtcTime.getMean());

	fprintf(fpSal, "Tiempo total: %.6f[s]\n", t2-t1);
	fprintf(fpSal, "Tiempo por cada par: %.6f[s]\n", (t2-t1)/(nEntren*(double)nPrueba));
	fprintf(fpSal, "Tiempo optimizado por cada par: %.6f[s]\n", rec.ti.calTime.getMean()+rec.ti.gtcTime.getMean());

	fclose(fpSal);
	// Destructores
	delete[] imEntren;
	delete[] imPrueba;
	return 0;
}


int L_ProgramasRecPuntInt::main_CalcTransf(int argc, char *argv[])
{
	if (argc != 5)
	{
		printf("calctr imagenRef.bmp imagenPru.bmp param.txt output.txt\n");
		return 1;
	}
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	L_KdTreeBBF tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	FILE *fp;
	L_ImageGrayDouble imRef, imPru;
	//Punteros para recorrer detecciones
	L_TransfPuntInt2DNodoPtr *ptrTr;
	L_CalceNodo *ptrCal;
	L_TransfAfinPI2D *ptrAf;

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.leeParametrosPropaga(argv[3]);
	rec.grabaParametrosPropaga(argv[3]);

	if ( imRef.readImage(argv[1]) == false )
	{
		printf("No se puede read la imagen de referencia %s\n", argv[1]);
		return 1;
	}
	if ( imPru.readImage(argv[2]) == false )
	{
		printf("No se puede read la imagen de test %s\n", argv[2]);
		return 1;
	}

	rec.agregaImagenEntrenRobaObjetoDe(imRef);
	rec.finalizaEntren();
	rec.procesaImagenesPrueba(imPru);

	fp = fopen(argv[4],"w");
	if (fp==NULL)
	{
		printf("No se puede escribir en archivo de output %s\n", argv[4]);
		return 1;
	}
	fprintf(fp,"%ld\n", long(rec.traL.size())); // Num transf
	for (ptrTr = rec.traL.root; ptrTr!=NULL; ptrTr=ptrTr->sig)
	{
		ptrAf=(L_TransfAfinPI2D *)(ptrTr->c);
		if (ptrAf == NULL)
		{
			printf("Problema interno de casting\n");
			return 1;
		}
		fprintf(fp, "%.4g %.4g %.4g\n", ptrAf->m11, ptrAf->m12, ptrAf->tx);
		fprintf(fp, "%.4g %.4g %.4g\n", ptrAf->m21, ptrAf->m22, ptrAf->ty);
		fprintf(fp, "%ld\n", long(ptrTr->c->calL.size()));
		for (ptrCal = ptrTr->c->calL.root; ptrCal != NULL ; ptrCal = ptrCal->sig)
		{
			fprintf(fp, "%.4g %.4g %.4g %.4g\n", ptrCal->c.dRef->x0, ptrCal->c.dRef->y0, ptrCal->c.dPru->x0, ptrCal->c.dPru->y0);
		}
	}
	fclose(fp);
	return 0;
}


///

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_leer_calibracion_camara(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_String str;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Parámetros cámara (*.parcam)", Fl_File_Chooser::SINGLE, "Leer parámetros cámara");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	str = dlg->value();
	if (obj->calib.leerParametrosCamara(str.c_str()) == false )
		printf("No se pudo read el archivo \"%s\"\n", str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_guardar_calibracion_camara(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_FileName nomarch;
	L_String str;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Parámetros cámara (*.parcam)", Fl_File_Chooser::CREATE, "Guardar parámetros cámara");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	nomarch.asign(dlg->value());
	if (nomarch.ext[0] == 0) // Sin extension
		nomarch.ext = "parcam";
	str = nomarch.dirnameext();
	if (obj->calib.grabarParametrosCamara(str.c_str()) == false )
		printf("No se pudo guardar el archivo \"%s\"\n", str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_leer_puntos(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_String str;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Puntos para calibración (*.puncam)", Fl_File_Chooser::SINGLE, "Leer puntos para calibración");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	str = dlg->value();
	if (obj->calib.leerArchivoPuntos(str.c_str()) == false )
		printf("No se pudo read el archivo \"%s\"\n", str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_guardar_puntos(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_FileName nomarch;
	L_String str;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Puntos para calibración (*.puncam)", Fl_File_Chooser::CREATE, "Guardar puntos para calibración");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	nomarch.asign(dlg->value());
	if (nomarch.ext[0] == 0) // Sin extension
		nomarch.ext = "puncam";
	str = nomarch.dirnameext();
	if (obj->calib.grabarArchivoPuntos(str.c_str()) == false )
		printf("No se pudo guardar el archivo \"%s\"\n", str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_salir(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->salir = true;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_resetear_calibracion_camara(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->calib.resetearCamara();
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_fovx_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double fovx=obj->calib.camInter.fovXCorregida();
	fovx = M_PI/180*L_PideNumFLTK("FovX[°]:", fovx*180/M_PI);
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), c->k1v(), c->k2v(), 0, 0, (c->resXv()/2)/tan(fovx/2), c->dfYv(), 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_fovy_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double fovy=obj->calib.camInter.fovYCorregida();
	fovy = M_PI/180*L_PideNumFLTK("FovY[°]:", fovy*180/M_PI);
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), c->k1v(), c->k2v(), 0, 0, c->dfXv(), (c->resYv()/2)/tan(fovy/2), 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_dfx_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double v;
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	v = L_PideNumFLTK("Distancia focal X[pix/tan]:", c->dfXv());
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), c->k1v(), c->k2v(), 0, 0, v, c->dfYv(), 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_dfy_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double v;
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	v = L_PideNumFLTK("Distancia focal Y[pix/tan]:", c->dfYv());
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), c->k1v(), c->k2v(), 0, 0, c->dfXv(), v, 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_k1_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double v;
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	v = L_PideNumFLTK("Constante q1[1/tan^2]:", c->k1v());
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), v, c->k2v(), 0, 0, c->dfXv(), c->dfYv(), 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_k2_manual(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	double v;
	L_CamaraDistorsionRadial *c = &obj->calib.camInter;
	v = L_PideNumFLTK("Constante q2[1/tan^4]:", c->k2v());
	c->fijaParametros(c->resXv(), c->resYv(), c->cenXv(), c->cenYv(), c->k1v(), v, 0, 0, c->dfXv(), c->dfYv(), 0);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_permitir_modif_k2(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->modificark2 = !obj->modificark2;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_referencia(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_String str;
	L_ImageRGBUchar imRef;
	L_ImageGrayDouble chica;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "(*.jpg)\t(*.jpeg)\t(*.bmp)", Fl_File_Chooser::SINGLE, "Leer imagen de referencia");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	str = dlg->value();
	if (obj->calib.fijaImagenReferencia(str.c_str()) == false )
		printf("No se pudo read el archivo \"%s\"\n", str.c_str());
	chica.cambiaTamanoPixelado(obj->calib.imRef, obj->venRef->vent.w(), obj->venRef->vent.h());
	imRef = chica;
	obj->venRef->dibuja(imRef);
	str.sprintf_cpy(1000, "Referencia %dx%d", obj->calib.imRef.lx, obj->calib.imRef.ly);
	obj->venRef->vent.label(str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_ancho_mm(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->calib.lxMM = L_PideNumFLTK("width imagen referencia[mm]:", obj->calib.lxMM);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_fijar_alto_mm(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->calib.lyMM = L_PideNumFLTK("Alto imagen referencia[mm]:", obj->calib.lxMM);
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_cambiar_camara(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->pedirCambioCamara = !obj->pedirCambioCamara;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_leer_config_referencia(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_String str;
	L_ImageRGBUchar imRef;
	L_ImageGrayDouble chica;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Configuracion de imagen referencia (*.cfgref)", Fl_File_Chooser::SINGLE, "Leer Configuracion de Referencia");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	str = dlg->value();
	if (obj->calib.leerArchivoConfiguracion(str.c_str()) == false )
		printf("No se pudo read el archivo \"%s\"\n", str.c_str());
	chica.cambiaTamanoPixelado(obj->calib.imRef, obj->venRef->vent.w(), obj->venRef->vent.h());
	imRef = chica;
	obj->venRef->dibuja(imRef);
	str.sprintf_cpy(1000, "Referencia %dx%d", obj->calib.imRef.lx, obj->calib.imRef.ly);
	obj->venRef->vent.label(str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_guardar_config_referencia(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	L_FileName nomarch;
	L_String str;
	Fl_File_Chooser *dlg = new Fl_File_Chooser(".", "Configuracion de imagen referencia (*.cfgref)", Fl_File_Chooser::CREATE, "Leer Configuracion de Referencia");
	dlg->show();
	while (dlg->shown()) // Esperando a que ingrese el archivo... .. .  .
		Fl::wait();
	if (dlg->value() == NULL) // Se presiono [cancel]
		return;
	nomarch.asign(dlg->value());
	if (nomarch.ext[0] == 0) // Sin extension
		nomarch.ext = "cfgref";
	str = nomarch.dirnameext();
	if (obj->calib.grabarArchivoConfiguracion(str.c_str()) == false )
		printf("No se pudo guardar el archivo \"%s\"\n", str.c_str());
	delete dlg;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_capturar_calces(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->pedirCapturaCalces = !obj->pedirCapturaCalces;
}

#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_borrar_calces(Fl_Widget *widget, void *param)
{
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->calib.puntos.resetear();
}

#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_CalibradorCamaraFLTK::cb_calibracion_camara_automatica(Fl_Widget *widget, void *param)
{
	// Aca la calibracion automatica usando TsaiCM...
	L_CalibradorCamaraFLTK *obj = (L_CalibradorCamaraFLTK *)param;
	obj->calib.calibrar_Tsai();
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
int L_CalibradorCamaraFLTK::main_calibrCam_obj(int argc, char *argv[])
{
	// Setear variables iniciales
	salir = false;
	estado = 0;
	modificark2 = false;
	pedirCapturaCalces=false;
	pedirCambioCamara = false;
	int cuadx;
	int cuady;

	Fl::visual(FL_RGB);
	Fl_Window *mainWin = new Fl_Window(0, 0, 320, 240, "Calibrador cámara");
	Fl_Menu_Item menuitems[] = {
		{"&Archivo", 0, 0, 0, FL_SUBMENU},
		{	"Leer calibración cámara", 0, &cb_leer_calibracion_camara, this},
		{	"Guardar calibración cámara", 0, &cb_guardar_calibracion_camara, this},
		{	"Leer puntos", 0, &cb_leer_puntos, this},
		{	"Guardar puntos", 0, &cb_guardar_puntos, this},
		{	"Salir", 0, &cb_salir, this},
		{	0},
		{"&Manual", 0, 0, 0, FL_SUBMENU},
		{	"Resetear calibración", 0, &cb_resetear_calibracion_camara, this},
		{	"Fijar FovX actual manual", 0, &cb_fijar_fovx_manual, this},
		{	"Fijar FovY actual manual", 0, &cb_fijar_fovy_manual, this},
		{	"Fijar dist focX manual", 0, &cb_fijar_dfx_manual, this},
		{	"Fijar dist focY manual", 0, &cb_fijar_dfy_manual, this},
		{	"Fijar k1 manual", 0, &cb_fijar_k1_manual, this},
		{	"Fijar k2 manual", 0, &cb_fijar_k2_manual, this},
		{	"Modificacion k2", 0, &cb_permitir_modif_k2, this},
		{	0},
		{"&Automática", 0, 0, 0, FL_SUBMENU},
		{   "Leer configuración de la referencia", 0, &cb_leer_config_referencia, this},
		{   "Guardar configuración de la referencia", 0, &cb_guardar_config_referencia, this},
		{   "Fijar imagen de referencia", 0, &cb_fijar_referencia, this},
		{   "Fijar width en papel", 0, &cb_fijar_ancho_mm, this},
		{   "Fijar height en papel", 0, &cb_fijar_alto_mm, this},
		{   "Cambiar cámara (estéreo)", 0, &cb_cambiar_camara, this},
		{	"Capturar calces", 0, &cb_capturar_calces, this},
		{	"clean calces", 0, &cb_borrar_calces, this},
		{	"Calibración automática", 0, &cb_calibracion_camara_automatica, this},
		{	0},
		{0}
	};
	Fl_Menu_Bar *menu = new Fl_Menu_Bar(0, 0, 320, 30);
	menu->copy(menuitems);
	Fl_Text_Buffer *buffer = new Fl_Text_Buffer(1000);
	Fl_Text_Display *editor = new Fl_Text_Display(0, 30, 600, 370);
	editor->buffer(buffer);
	mainWin->show();

	venCam = new L_VentanaImagenFLTK(400, 0, 320, 240, "Cámara");
	venRef = new L_VentanaImagenFLTK(0, 300, 320, 240, "Referencia");
	venCal = new L_VentanaImagenFLTK(400, 300, 320, 240, "Calces");

	L_ImageGrayDouble imx;
	L_ImageRGBUchar imTr;
	L_String mensaje;
	double dfMod, maxk1, maxk2;
	int i, j;

	L_CapturadorImagen capt;

	if (capt.capturarImagen() == false)
	{
		printf("No se puede capturar\n");
		return 1;
	}
	if (capt.im.lx == 0)
		printf("Error\n");

	if (capt.estereo.esEstereo)
		capt.estereo.usarIzq0Der1Ambas2 = 2;

	while(salir == false)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No se puede capturar\n");
			return 1;
		}
		if (pedirCambioCamara && capt.estereo.esEstereo)
			imx=capt.imDer;
		else
			imx=capt.im;

		calib.camInter.corrigeImagenBilineal(capt.im, imCorr);
		//calib.camInter.muestraImagenError(capt.im, imCorr, 4);

		// Dibujar cuadradito que muestra la distorsion como click
		dfMod = sqrt(calib.camInter.dfXv()*calib.camInter.dfXv() + calib.camInter.dfYv()*calib.camInter.dfYv());
		maxk1 = (10.0/500/500) * dfMod * dfMod;
		maxk2 = (10.0/500/500/500/500) * dfMod * dfMod * dfMod * dfMod;
		cuadx = (int)(  0.5*(2.0*calib.camInter.k1v() + maxk1)*imCorr.lx / maxk1  );
		cuady = (int)(  0.5*(2.0*calib.camInter.k2v() + maxk2)*imCorr.ly / maxk2  );
		for (j=cuady-1; j<=cuady+1; j++)
		{
			for (i=cuadx-1; i<=cuadx+1; i++)
			{
				if (i>=0 && i<imCorr.lx && j>=0 && j<imCorr.ly)
				{
					imCorr.pix(i,j,0)=255;
					imCorr.pix(i,j,1)=255;
					imCorr.pix(i,j,2)=255;
				}
			}
		}
		venCam->dibuja(imCorr); // Mostrar la corregida; asi hay realimentacion al calibrar manualmente

		switch(estado)
		{
		case 0: // Descanso
			if (venCam->mouseE()!=L_MouseLeido)
			{
				double k1, k2;
				venCam->mouseE()=L_MouseLeido;
				k1 = maxk1 * (venCam->mouseX()-imCorr.lx/2) / (double) imCorr.lx;
				k2 = maxk2 * (venCam->mouseY()-imCorr.ly/2) / (double) imCorr.ly;
				if (modificark2)
					calib.camInter.fijaDistorsionK1_Tsai(k1);
				else
					calib.camInter.fijaDistorsionK1_Tsai(k1);
				estado = 2;
			}
			if (pedirCambioCamara)
			{
				pedirCambioCamara = false;
				if (capt.estereo.esEstereo)
					capt.estereo.usarIzq0Der1Ambas2 = 1-capt.estereo.usarIzq0Der1Ambas2;
			}
			if (pedirCapturaCalces == true)
			{
				pedirCapturaCalces = false;
				estado = 1;
			}
			break;
		case 1: // Calibracion automatica
			if (venCal->mouseE()!=L_MouseLeido)
			{
				venCal->mouseE()=L_MouseLeido;
				estado = 0;
			}
			if (calib.imRef.data() == NULL)
			{
				fl_message("No se eligio una imagen de referencia");
				estado = 0;
				break;
			}
			if (calib.rpiEntrenado == false)
				calib.procesaImagenReferencia();
			calib.procesaNuevaImagen(imx);
			if (calib.rpi.traL.size() > 0)
				mensaje.sprintf_cpy(1000, "Calces: %d de %d\n", calib.rpi.traL.root->c->calL.size(), calib.rpi.calL.size());
			else
				mensaje.sprintf_cpy(1000, "Calces: 0 de %d\n", calib.rpi.calL.size());
			calib.rpi.dibujaTransformaciones(imTr, 1, 1, 1, 1);
			//calib.rpi.calL.dibuja(imTr);
			//venCal->label(mensaje.elem);
			venCal->resize(imTr.lx, imTr.ly);
			venCal->dibuja(imTr);
			break;

		case 2: // Calibracion manual
			if (venCam->mouseE()!=L_MouseLeido)
			{
				venCam->mouseE()=L_MouseLeido;
				double k1, k2;
				k1 = maxk1 * (venCam->mouseX()-imCorr.lx/2) / (double) imCorr.lx;
				k2 = maxk2 * (venCam->mouseY()-imCorr.ly/2) / (double) imCorr.ly;
				if (modificark2)
					calib.camInter.fijaDistorsionK1_Tsai(k1);
				else
					calib.camInter.fijaDistorsionK1_Tsai(k1);
				estado = 2;
			}
			if (pedirCambioCamara)
			{
				pedirCambioCamara = false;
				if (capt.estereo.esEstereo)
					capt.estereo.usarIzq0Der1Ambas2 = 1-capt.estereo.usarIzq0Der1Ambas2;
			}
			if (pedirCapturaCalces == true)
			{
				pedirCapturaCalces = false;
				//calib.camInter.fijaDistorsionQ1Q2(calib.camInter.q1, calib.camInter.q2);
				estado = 1;
			}
			break;

		default:
			printf("Error interno\n");
			return 1;
		}

		mensaje = "";
		mensaje.sprintf_cat(1000, "resX = %d[pix]   cenX = %.2f[pix]\n", calib.camInter.resXv(), calib.camInter.cenXv());
		mensaje.sprintf_cat(1000, "resY = %d[pix]   cenY = %.2f[pix]\n", calib.camInter.resYv(), calib.camInter.cenYv());
		mensaje.sprintf_cat(1000, "dfX  = %.3f[pix/tan]   fovX=%.2f[°]\n", calib.camInter.dfXv(), calib.camInter.fovXCorregida()*180/M_PI);
		mensaje.sprintf_cat(1000, "dfY  = %.3f[pix/tan]   fovY=%.2f[°]\n", calib.camInter.dfYv(), calib.camInter.fovYCorregida()*180/M_PI);
		mensaje.sprintf_cat(1000, "k1   = %.5f[1/tan^2]\n", calib.camInter.k1v());
		mensaje.sprintf_cat(1000, "k2   = %.5f[1/tan^4]\n", calib.camInter.k2v());
		if (modificark2)
			mensaje.sprintf_cat(1000, "modificar k2 activado\n");
		else
			mensaje.sprintf_cat(1000, "modificar k2 desactivado\n");
		mensaje.sprintf_cat(1000, "-------------------\n");
		mensaje.sprintf_cat(1000, "imRef: %d[pix] x %d[pix]\n", calib.imRef.lx, calib.imRef.ly);
		mensaje.sprintf_cat(1000, "imRef: %.1f[mm] x %.1f[mm]\n", calib.lxMM, calib.lyMM);
		mensaje.sprintf_cat(1000, "nCuadros = %d\n", calib.puntos.size());
		mensaje.sprintf_cat(1000, "(%.2f calces mean)\n", calib.puntos.nPuntosEstad.getMean());
		buffer->text(mensaje.c_str());
		Fl::check();
	}
	delete mainWin;
	delete venCam;
	delete venRef;
	delete venCal;
	return 0;
}
#endif // __COMPAT_FLTK__


#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_captStereo(int argc, char *argv[])
{
	L_CapturadorImagen capt;
	L_VentanaImagen venIzq(0, 0, 320, 240, "Izq");
	L_VentanaImagen venDer(400, 0, 320, 240, "Der");
	L_VentanaImagen venIzq2(0, 300, 320, 240, "Izq2");
	L_VentanaImagen venDer2(400, 300, 320, 240, "Der2");
	L_ImageRGBUchar imL, imR;
	//L_Camara_RadialTang modIzq, modDer;
	L_CamaraDistorsionRadial modIzq, modDer;
	modIzq.fijaParametrosVidere_lentesOrig(); //.fijaParametrosVidereIzq();
	modDer.fijaParametrosVidere_lentesOrig(); //.fijaParametrosVidereDer();
	char str[100];
	int i = 0;
	if (capt.capturarImagen() == false)
	{
		printf("No se puede capturar\n");
		return 1;
	}
	if (capt.estereo.esEstereo == false)
	{
		printf("Se requiere camara estereo\n");
		return 2;
	}
	capt.estereo.usarIzq0Der1Ambas2 = 2;
	while (capt.capturarImagen())
	{
		venIzq.dibuja(capt.im);
		venDer.dibuja(capt.imDer);
		modIzq.corrigeImagenBilineal(capt.im, imL);
		modDer.corrigeImagenBilineal(capt.imDer, imR);
		venIzq2.dibuja(imL);
		venDer2.dibuja(imR);
		venIzq.check();
		switch(venIzq.mouseE())
		{
		case L_MousePresiona:
			venIzq.mouseE() = L_MouseLeido;
			sprintf(str,"left%d.bmp", i);
			capt.im.saveImage(str);
			sprintf(str,"right%d.bmp", i);
			capt.imDer.saveImage(str);
			printf("%s\n", str);
			i++;
			break; // Error: faltaba esto
		default:
			venIzq.mouseE()=L_MouseLeido;
		};
	}
	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)



#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camRecPuntInt(int argc, char *argv[])
{
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	//L_PDoG_SURFG sift;
	//L_PDoG_SURFG sift; sift.genPiram.usasqrt2 = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_SURF sift;
	//L_PHarris_SIFT sift;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = true; sift.sgen.orig = true; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = false; sift.sgen.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_OpenSURF sift;

	//const char *nomArchParams = "paramsRecPuntCam_OpenSURF.txt";
	const char *nomArchParams = "paramsRecPuntCam_SURF_Lonco.txt";
	//const char *nomArchParams = "paramsRecPuntCam_SIFT.txt";

	L_KdTreeBBF tree; tree.nCompMax = 64; tree.usePercentage = false;
	//L_Flann tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;
	L_CapturadorImagen capt;
	L_DescriptorLista desLRef;
	L_ShapeArray lins, linsCamara;
	int x0=0, y0=0;
	bool entrenado = false;
	L_Array<const L_ImageGrayDouble *> imarr;

	//capt.fijarRedimensionar(true, 320, 240);

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	int X = im.lx + 100;
	int Y = im.ly + 100;
	int Z = 2*im.ly + 200;

#ifdef __COMPAT_FLTK__
	int anchoPant = Fl::w();
	int altoPant = Fl::h();
#else __COMPAT_FLTK__
	int anchoPant = 800;
	int altoPant = 600;
#endif

	if (X + im.lx > anchoPant)
		X = anchoPant - im.lx;

	if (Y + im.ly > altoPant)
		Y = altoPant - im.ly - 50;

	if (Z + im.ly > altoPant)
		Z = altoPant - im.ly - 10;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venPunRef(X, 0, im.lx, im.ly, "Puntos ref (click=salir)");
	L_VentanaImagen venPunPru(0, Y, im.lx, im.ly, "Puntos test");
	L_VentanaImagen venCal(X, Y, im.lx, im.ly, "Calces");
	L_VentanaImagen venCalRec(0, Z, im.lx, im.ly, "Calces Mejor Transformacion");
	L_VentanaImagen venTrRec(X, Z, im.lx, im.ly, "Transformaciones");

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.leeParametrosPropaga(nomArchParams);
	rec.grabaParametrosPropaga(nomArchParams);
	rec.debugActivo(true);

	while(true)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);

		if (venCam.mouseE() == L_MousePresiona)
		{
			venCam.mouseE() = L_MouseLeido;
			x0 = venCam.mouseX();
			y0 = venCam.mouseY();
			entrenado = false; // Para subir el frame rate
		}
		if (venCam.mouseE() == L_MouseArrastra)
		{
			venCam.mouseE() = L_MouseLeido;
			entrenado = false; // Para subir el frame rate
			linsCamara.drawRectangle(x0, y0, venCam.mouseX(), venCam.mouseY());
		}
		if (venCam.mouseE() == L_MouseSuelta)
		{
			venCam.mouseE() = L_MouseLeido;
			if (venCam.mouseX()-x0 > 0 && venCam.mouseY()-y0 > 0) // Por si hay problemas de sincronismo
			{
				int xM=venCam.mouseX()-x0+1, yM=venCam.mouseY()-y0+1;
				if (xM > capt.im.lx - x0)
					xM = capt.im.lx - x0 -1;
				if (yM > capt.im.ly - y0)
					yM = capt.im.ly - y0 -1;
				imTrozo.subImageOf(capt.im, x0, y0, xM, yM);
				imTrozox = imTrozo;
				rec.reseteaEntren();
				rec.agregaImagenEntrenCopiaObjetoDe(imTrozox);
				desLRef.clear();
				rec.desLRef.copyListOn(desLRef);
				rec.finalizaEntren();
				entrenado = true;
				X = im.lx+imTrozo.lx + 100;

				if (X + im.lx + imTrozo.lx > anchoPant)
					X = anchoPant - im.lx - imTrozo.lx;

				venCal.resize(X, Y, im.lx+imTrozo.lx, im.ly);
				venCalRec.resize(0, Z, im.lx+imTrozo.lx, im.ly);
				venTrRec.resize(X, Z, im.lx+imTrozo.lx, im.ly);
			}
		}
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibuja(imTmp);
		if (entrenado)
		{
			rec.procesaImagenesPrueba(imx);

			printf("%d descriptores\n", rec.desLPru.size());

			rec.generaArregloImagenesParaCalces(imarr);

			desLRef.dibujaFlechas(lins);
			imTmp = imTrozox;
			imTmp.genDrawing(lins);
			lins.resize(0);
			venPunRef.dibuja(imTmp);

			rec.desLPru.dibujaFlechas(lins);
			imTmp = im;
			imTmp.genDrawing(lins);
			lins.resize(0);
			venPunPru.dibuja(imTmp);

			rec.calL.genLineas(lins);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venCal.dibuja(imTmp);

			rec.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venCalRec.dibuja(imTmp);
			rec.traL.dibujaTransformaciones(lins, 0, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venTrRec.dibuja(imTmp);

			rec.ti.escribeResumen(stdout);
			rec.reseteaPrueba();
		}

		if (venCal.mouseE() != L_MouseLeido || venPunRef.mouseE() != L_MouseLeido  || venPunPru.mouseE() != L_MouseLeido  || venCalRec.mouseE() != L_MouseLeido  || venTrRec.mouseE() != L_MouseLeido)
			return 0;
		venCal.check();
	}
	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)



#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camRecPuntInt_pose_plana(int argc, char *argv[])
{
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	//L_PDoG_SURFG sift;
	//L_PDoG_SURFG sift; sift.genPiram.usasqrt2 = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_SURF sift;
	//L_PHarris_SIFT sift;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = true; sift.sgen.orig = true; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = false; sift.sgen.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_OpenSURF sift;

	L_KdTreeBBF tree; tree.nCompMax = 64; tree.usePercentage = false;
	//L_Flann tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	L_Pose3D_cuat pose;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;
	L_CapturadorImagen capt;
	L_DescriptorLista desLRef;
	L_ShapeArray lins, linsCamara;
	int x0=0, y0=0;
	bool entrenado = false;
	L_Array<const L_ImageGrayDouble *> imarr;

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	int X = im.lx + 100;
	int Y = im.ly + 100;
	int Z = 2*im.ly + 200;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venPunRef(X, 0, im.lx, im.ly, "Puntos ref (click=salir)");
	L_VentanaImagen venPunPru(0, Y, im.lx, im.ly, "Puntos test");
	L_VentanaImagen venCal(X, Y, im.lx, im.ly, "Calces");
	L_VentanaImagen venCalRec(0, Z, im.lx, im.ly, "Calces Mejor Transformacion");
	L_VentanaImagen venTrRec(X, Z, im.lx, im.ly, "Transformaciones");

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.debugActivo(true);

	while(true)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);

		if (venCam.mouseE() == L_MousePresiona)
		{
			venCam.mouseE() = L_MouseLeido;
			x0 = venCam.mouseX();
			y0 = venCam.mouseY();
			entrenado = false; // Para subir el frame rate
		}
		if (venCam.mouseE() == L_MouseArrastra)
		{
			venCam.mouseE() = L_MouseLeido;
			entrenado = false; // Para subir el frame rate
			linsCamara.drawRectangle(x0, y0, venCam.mouseX(), venCam.mouseY());
		}
		if (venCam.mouseE() == L_MouseSuelta)
		{
			venCam.mouseE() = L_MouseLeido;
			if (venCam.mouseX()-x0 > 0 && venCam.mouseY()-y0 > 0) // Por si hay problemas de sincronismo
			{
				int xM=venCam.mouseX()-x0+1, yM=venCam.mouseY()-y0+1;
				if (xM > capt.im.lx - x0)
					xM = capt.im.lx - x0 -1;
				if (yM > capt.im.ly - y0)
					yM = capt.im.ly - y0 -1;
				imTrozo.subImageOf(capt.im, x0, y0, xM, yM);
				imTrozox = imTrozo;
				rec.reseteaEntren();
				rec.agregaImagenEntrenCopiaObjetoDe(imTrozox);
				desLRef.clear();
				rec.desLRef.copyListOn(desLRef);
				rec.finalizaEntren();
				entrenado = true;
				X = im.lx+imTrozo.lx + 100;
				venCal.resize(X, Y, im.lx+imTrozo.lx, im.ly);
				venCalRec.resize(0, Z, im.lx+imTrozo.lx, im.ly);
				venTrRec.resize(X, Z, im.lx+imTrozo.lx, im.ly);
			}
		}
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibuja(imTmp);
		if (entrenado)
		{
			rec.procesaImagenesPrueba(imx);

			printf("%d descriptores\n", rec.desLPru.size());

			rec.generaArregloImagenesParaCalces(imarr);

			desLRef.dibujaFlechas(lins);
			imTmp = imTrozox;
			imTmp.genDrawing(lins);
			lins.resize(0);
			venPunRef.dibuja(imTmp);

			rec.desLPru.dibujaFlechas(lins);
			imTmp = im;
			imTmp.genDrawing(lins);
			lins.resize(0);
			venPunPru.dibuja(imTmp);

			rec.calL.genLineas(lins);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venCal.dibuja(imTmp);

			rec.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venCalRec.dibuja(imTmp);
			rec.traL.dibujaTransformaciones(lins, 0, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			venTrRec.dibuja(imTmp);

			rec.ti.escribeResumen(stdout);
			rec.reseteaPrueba();
		}

		if (venCal.mouseE() != L_MouseLeido || venPunRef.mouseE() != L_MouseLeido  || venPunPru.mouseE() != L_MouseLeido  || venCalRec.mouseE() != L_MouseLeido  || venTrRec.mouseE() != L_MouseLeido)
			return 0;
		venCal.check();
	}
	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)






#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camRecPuntInt_video_0(int argc, char *argv[])
{
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	//L_PDoG_SURFG sift;
	//L_PDoG_SURFG sift; sift.genPiram.usasqrt2 = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_SURF sift;
	//L_PHarris_SIFT sift;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = true; sift.sgen.orig = true; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = false; sift.sgen.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_OpenSURF sift;

	L_KdTreeBBF tree; //tree.usePercentage = false;
	//L_BusqSecDescr tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;

	std::vector<double> xp, yp;
	FILE *fp = NULL;


	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\heman.avi");
	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\casteloEspecial.avi");
	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\invitation.avi");
	L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\nti.avi");

	const char name[] = "matches_nti_tr.txt";

	L_DescriptorLista desLRef;
	L_ShapeArray lins, linsCamara;
	int x0=0, y0=0;
	bool entrenado = false;
	L_Array<const L_ImageGrayDouble *> imarr;

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	int X = im.lx + 100;
	int Y = 300;//im.ly + 100;
	int Z = 600;//2*im.ly + 200;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venPunRef(X, 0, im.lx, im.ly, "Puntos ref (click=salir)");
	L_VentanaImagen venPunPru(0, Y, im.lx, im.ly, "Puntos test");
	L_VentanaImagen venCal(X, Y, im.lx, im.ly, "Calces");
	L_VentanaImagen venCalRec(0, Z, im.lx, im.ly, "Calces Mejor Transformacion");
	L_VentanaImagen venTrRec(X, Z, im.lx, im.ly, "Transformaciones");

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.debugActivo(true);
	rec.ti.escribeResumen(stdout);

	while(true)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibujaRedimensiona(imTmp);


		// Formato:  x der,  y aba
		//
		//  numFeatures numFrames
		//  x  y  1    x  y  1    x  y  1     -> num frames
		// -1 -1 -1    x  y  1    x  y  1
		//  |
		//  V  num features

		// Esta funcion write los datos al reves

		if (!entrenado)
		{
			imTrozox = capt.im;
			rec.reseteaEntren();
			rec.agregaImagenEntrenCopiaObjetoDe(imTrozox);
			xp.resize(rec.desLRef.size());
			yp.resize(rec.desLRef.size());
			rec.desLRef.fijaIdCorrelativos();

			desLRef.clear();
			rec.desLRef.copyListOn(desLRef);
			rec.finalizaEntren();
			entrenado = true;
			X = im.lx+imTrozo.lx + 100;
			venCal.resize(X, Y, im.lx+imTrozo.lx, im.ly);
			venCalRec.resize(0, Z, im.lx+imTrozo.lx, im.ly);
			venTrRec.resize(X, Z, im.lx+imTrozo.lx, im.ly);
			fp = fopen(name, "w");
			continue;
		}
		else
		{
			rec.procesaImagenesPrueba(imx);
			rec.ti.escribeResumen(stdout);
			printf("%d calces test\n", rec.desLPru.size());

			rec.generaArregloImagenesParaCalces(imarr);

			desLRef.dibujaFlechas(lins);
			imTmp = imTrozox;
			imTmp.genDrawing(lins);
			lins.resize(0);
			imTmp.halfResolution();
			venPunRef.dibujaRedimensiona(imTmp);

			rec.desLPru.dibujaFlechas(lins);
			imTmp = im;
			imTmp.genDrawing(lins);
			lins.resize(0);
			imTmp.halfResolution();
			venPunPru.dibujaRedimensiona(imTmp);

			rec.calL.genLineas(lins);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venCal.dibujaRedimensiona(imTmp);

			rec.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venCalRec.dibuja(imTmp);

			rec.traL.dibujaTransformaciones(lins, 0, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venTrRec.dibujaRedimensiona(imTmp);

			// Guardar datos

			if (rec.traL.size() != 0)
			{
				L_CalceNodo *cal, *root;
				for (int i=0; i<(int)xp.size(); i++)
					{xp[i]=0.0; yp[i]=0.0;}
				root = rec.traL.root->c->calL.root;
				for (cal = root; cal!=NULL; cal=cal->sig)
				{
					xp[cal->c.dRef->id] = cal->c.dPru->x0;
					yp[cal->c.dRef->id] = cal->c.dPru->y0;
				}
				for (int i=0; i<(int)xp.size(); i++)
				{
					if (xp[i] != 0)
						fprintf(fp, "%.3f %.3f 1 ", xp[i], yp[i]);  // x=der, y=aba, (0,0) = esquina
					else
						fprintf(fp, "-1 -1 -1 ");  // top left
				}
				fprintf(fp, "\n");
			}
			rec.reseteaPrueba();
		}

		if (venCal.mouseE() != L_MouseLeido || venPunRef.mouseE() != L_MouseLeido  || venPunPru.mouseE() != L_MouseLeido  || venCalRec.mouseE() != L_MouseLeido  || venTrRec.mouseE() != L_MouseLeido)
			return 0;
		venCal.check();
	}

	if (fp != NULL)
		{fclose(fp); fp = NULL;}


	//L_Matrix m, m2;
	//m.readFile("C:\\Users\\ploncomi\\Documents\\matches_invitation_tr.txt"); // m: aba = nfr, der = 3*nfeat
	//m2.reallocate(m.lj/3, m.li*3);  // m2: aba = nfeat, der = 3*nfr
	//for (int i=0; i<m.li; i++)  // i: nFr
	//	for (int j=0; j<m2.li; j++) // j: nFeat
	//		for (int k=0; k<3; k++)
	//			m2[j][3*i+k] = m[i][3*j+k];
	//FILE *fp = fopen("C:\\Users\\ploncomi\\Documents\\matches_invitation.txt", "w");
	//if (fp == NULL)
	//	return 1;
	//fprintf(fp, "%d %d\n", m.li, m.lj/3);
	//m2.saveFile(fp);
	//fclose(fp);


	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)


#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camRecPuntInt_video_1(int argc, char *argv[])
{
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	//L_PDoG_SURFG sift;
	//L_PDoG_SURFG sift; sift.genPiram.usasqrt2 = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_SURF sift;
	//L_PHarris_SIFT sift;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = true; sift.sgen.orig = true; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = true; sift.sdirL.orig = false; sift.sgen.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = false;
	//L_PDoG_SIFT sift; sift.genPiram.usasqrt2 = false; sift.sgen.orig = false; sift.sdirL.orig = false; sift.genPiram.dupl = true; sift.val_minP_DoG = 0.01;
	//L_OpenSURF sift;

	L_KdTreeBBF tree; //tree.usePercentage = false;
	//L_BusqSecDescr tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;

	std::vector<double> xp, yp;
	L_Array<std::vector<double> > m;


	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\heman.avi");
	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\casteloEspecial.avi");
	//L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\invitation.avi");
	L_CapturadorImagen capt("C:\\Users\\ploncomi\\Documents\\nti.avi");

	const char name[] = "C:\\Users\\ploncomi\\Documents\\matches_nti.txt";

	L_DescriptorLista desLRef;
	L_ShapeArray lins, linsCamara;
	int x0=0, y0=0;
	bool entrenado = false;
	L_Array<const L_ImageGrayDouble *> imarr;

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	int X = im.lx + 100;
	int Y = 300;//im.ly + 100;
	int Z = 600;//2*im.ly + 200;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venPunRef(X, 0, im.lx, im.ly, "Puntos ref (click=salir)");
	L_VentanaImagen venPunPru(0, Y, im.lx, im.ly, "Puntos test");
	L_VentanaImagen venCal(X, Y, im.lx, im.ly, "Calces");
	L_VentanaImagen venCalRec(0, Z, im.lx, im.ly, "Calces Mejor Transformacion");
	L_VentanaImagen venTrRec(X, Z, im.lx, im.ly, "Transformaciones");

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.debugActivo(true);
	rec.ti.escribeResumen(stdout);

	while(true)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			break;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibujaRedimensiona(imTmp);


		// Formato:  x der,  y aba
		//
		//  numFeatures numFrames
		//  x  y  1    x  y  1    x  y  1     -> num frames
		// -1 -1 -1    x  y  1    x  y  1
		//  |
		//  V  num features

		// Esta funcion guarda los datos al reves inicialmente

		if (!entrenado)
		{
			imTrozox = capt.im;
			rec.reseteaEntren();
			rec.agregaImagenEntrenCopiaObjetoDe(imTrozox);
			xp.resize(rec.desLRef.size());
			yp.resize(rec.desLRef.size());
			rec.desLRef.fijaIdCorrelativos();

			for (int i=0; i<(int)xp.size(); i++)
				xp[i] = -1;

			for (L_DescriptorNodo *des = rec.desLRef.root; des!=NULL; des=des->sig)
			{
				xp[des->c.id] = des->c.x0;
				yp[des->c.id] = des->c.y0;
			}

			// m: m[frame][3*nfeat]
			m.resize_swapping(m.size()+1);
			m[m.size()-1].resize(3*rec.desLRef.size());

			for (int i=0; i<(int)xp.size(); i++)
			{
				if (xp[i] != -1)
				{
					m[m.size()-1][3*i+0] = xp[i];
					m[m.size()-1][3*i+1] = yp[i];
					m[m.size()-1][3*i+2] = 1;
				}
				else
				{
					m[m.size()-1][3*i+0] = -1;
					m[m.size()-1][3*i+1] = -1;
					m[m.size()-1][3*i+2] = -1;
				}
			}

			desLRef.clear();
			rec.desLRef.copyListOn(desLRef);
			rec.finalizaEntren();
			entrenado = true;

			X = im.lx+imTrozo.lx + 100;
			venCal.resize(X, Y, im.lx+imTrozo.lx, im.ly);
			venCalRec.resize(0, Z, im.lx+imTrozo.lx, im.ly);
			venTrRec.resize(X, Z, im.lx+imTrozo.lx, im.ly);
			continue;
		}
		else
		{
			rec.procesaImagenesPrueba(imx);
			rec.ti.escribeResumen(stdout);
			printf("%d calces test\n", rec.desLPru.size());

			rec.generaArregloImagenesParaCalces(imarr);

			desLRef.dibujaFlechas(lins);
			imTmp = imTrozox;
			imTmp.genDrawing(lins);
			lins.resize(0);
			imTmp.halfResolution();
			venPunRef.dibujaRedimensiona(imTmp);

			rec.desLPru.dibujaFlechas(lins);
			imTmp = im;
			imTmp.genDrawing(lins);
			lins.resize(0);
			imTmp.halfResolution();
			venPunPru.dibujaRedimensiona(imTmp);

			rec.calL.genLineas(lins);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venCal.dibujaRedimensiona(imTmp);

			rec.traL.dibujaTransformaciones(lins, 1, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venCalRec.dibuja(imTmp);

			rec.traL.dibujaTransformaciones(lins, 0, 1, 1, 1);
			imTmp = im;
			imTmp.genDrawingMatches(imarr, lins);
			lins.resize(0);
			imTmp.halfResolution();
			venTrRec.dibujaRedimensiona(imTmp);

			// Guardar datos

			if (rec.traL.size() != 0)
			{
				for (int i=0; i<(int)xp.size(); i++)
					xp[i] = -1;

				L_CalceNodo *root = rec.traL.root->c->calL.root;
				for (L_CalceNodo *cal = root; cal!=NULL; cal=cal->sig)
				{
					xp[cal->c.dRef->id] = cal->c.dPru->x0;
					yp[cal->c.dRef->id] = cal->c.dPru->y0;
				}
				// m: m[frame][3*nfeat]
				m.resize_swapping(m.size()+1);
				m[m.size()-1].resize(3*xp.size());

				for (int i=0; i<(int)xp.size(); i++)
				{
					if (xp[i] != -1)
					{
						m[m.size()-1][3*i+0] = xp[i];
						m[m.size()-1][3*i+1] = yp[i];
						m[m.size()-1][3*i+2] = 1;
					}
					else
					{
						m[m.size()-1][3*i+0] = -1;
						m[m.size()-1][3*i+1] = -1;
						m[m.size()-1][3*i+2] = -1;
					}
				}
			}
			rec.reseteaPrueba();
		}

		if (venCal.mouseE() != L_MouseLeido || venPunRef.mouseE() != L_MouseLeido  || venPunPru.mouseE() != L_MouseLeido  || venCalRec.mouseE() != L_MouseLeido  || venTrRec.mouseE() != L_MouseLeido)
			return 0;
		venCal.check();
	}


	L_Matrix mat;

	// m:   m[frame][3*nfeat]
	// mat: mat[nfeat][3*nframe]

	mat.reallocate((int)m[0].size()/3, (int)m.size()*3);

	for (int i=0; i<mat.lj/3; i++)  // i: nFr
		for (int j=0; j<mat.li; j++) // j: nFeat
			for (int k=0; k<3; k++)
				mat(j,3*i+k) = m[i][3*j+k];

	FILE *fp = NULL;
	fp = fopen(name, "w");

	if (fp == NULL)
	{
		printf("No se pudo abrir %s para escritura\n", name);
		return 1;
	}

	fprintf(fp, "%d %d\n", mat.li, mat.lj/3);  // nFeatures  nFrames
	mat.saveFile(fp);  // mat[nfeat][3*nframe]

	fclose(fp);


	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)



#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_trackCaja(int argc, char *argv[])
{
#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	//L_OpenSURF sift;

	L_KdTreeBBF tree; tree.nCompMax = 64; tree.usePercentage = false;
	//L_Flann tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;
	L_CapturadorImagen capt("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Bases de Datos\\Brasil\\Tracking3D\\KBRealSequence.avi");
	L_DescriptorLista desLRef;
	L_ShapeArray lins;
	L_Descriptor3DArreglo desArr;
	L_CamaraPinhole cam;
	L_Pose3D_cuat pose;
	L_Pose3D ptr;
	L_CoordsCart3D pos, normal, horz;
	L_DescriptorNodo *des;
	L_CalceNodo *cal;
	std::vector<L_RayoCamara> arrRay;
	std::vector<L_CoordsCart3D> arrPos;
	L_RayoCamara ray;
	L_LevenbergMarquardt opt;
	double errAng = 45*M_PI/180;
	double distCam;
	int i;
	int x0[4], y0[4];

	for (i=0; i<500; i++)
		capt.capturarImagen();

	L_ImageGrayDouble ims[4];

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	double asp = 0.8;

	ims[0].readImage("lado0.png");
	ims[1].readImage("lado1.png");
	ims[2].readImage("lado2.png");
	ims[3].readImage("lado3.png");

	cam.fijaParametrosPhillips();

	// Entrenamiento caja
	// Calcular distancia a la caja
	std::vector<L_RayoCamara> r4p(4);
	std::vector<L_CoordsCart3D> p3p(4);
	cam.calcRayo(r4p[0], 0, 0);
	cam.calcRayo(r4p[1], 147*asp, 0);
	cam.calcRayo(r4p[2], 147*asp, 167);
	cam.calcRayo(r4p[3], 147*asp, 0);
	p3p[0].fijaXYZ(0, 0, 0);
	p3p[1].fijaXYZ(0, -147*asp, 0);
	p3p[2].fijaXYZ(0, -147*asp, -167);
	p3p[3].fijaXYZ(0, -0, -167);

	pose.calcPose_movido_a_fijo_4rayos(p3p, r4p);
	distCam = fabs(pose.pos.x); // Deberia quedar con signo negativo

	desArr.reserve(6000);
	// Agregar cara de adelante
	sift.calcDescrPuntInt(ims[0]);
	normal.x = -1;
	normal.y = 0;
	normal.z = 0;
	horz.x = 0;
	horz.y = -1;
	horz.z = 0;
	for (des = sift.desL.root; des!=NULL; des=des->sig)
	{
		//pos.x = 0.05 * (167 - des->c.y0) * asp; // 0*asp
		pos.x = 0;
		pos.y = -des->c.x0*asp;
		pos.z = -des->c.y0;
		desArr.resize(desArr.size()+1);
		desArr[desArr.size()-1].crear(des->c, pos, normal, horz, distCam, cam);
	}
	sift.desL.clear();

	// Agregar cara del lado derecho
	sift.calcDescrPuntInt(ims[1]);
	normal.x = 0;
	normal.y = -1;
	normal.z = 0;
	horz.x = 1;
	horz.y = 0;
	horz.z = 0;
	for (des = sift.desL.root; des!=NULL; des=des->sig)
	{
		pos.x = des->c.x0*asp;
		pos.y = -147.0*asp;
		pos.z = -des->c.y0;
		desArr.resize(desArr.size()+1);
		desArr[desArr.size()-1].crear(des->c, pos, normal, horz, distCam, cam);
	}
	sift.desL.clear();

	// Agregar cara de atras
	sift.calcDescrPuntInt(ims[2]);
	normal.x = 1;
	normal.y = 0;
	normal.z = 0;
	horz.x = 0;
	horz.y = 1;
	horz.z = 0;
	for (des = sift.desL.root; des!=NULL; des=des->sig)
	{
		//pos.x = 85*asp;
		pos.x = 85*asp - 0.05*des->c.y0;
		pos.y = (-147.0 + des->c.x0)*asp;
		pos.z = -des->c.y0;
		desArr.resize(desArr.size()+1);
		desArr[desArr.size()-1].crear(des->c, pos, normal, horz, distCam, cam);
	}
	sift.desL.clear();

	// Agregar cara del lado izquierdo
	sift.calcDescrPuntInt(ims[3]);
	normal.x = 0;
	normal.y = 1;
	normal.z = 0;
	horz.x = -1;
	horz.y = 0;
	horz.z = 0;
	for (des = sift.desL.root; des!=NULL; des=des->sig)
	{
		pos.x = (85 - des->c.x0)*asp;
		pos.y = 0*asp;
		pos.z = -des->c.y0;
		desArr.resize(desArr.size()+1);
		desArr[desArr.size()-1].crear(des->c, pos, normal, horz, distCam, cam);
	}
	sift.desL.clear();

	// Definir pose inicial
	pose.fijaCero();
	pose.pos.x = distCam; // distCam o 200; // Una estimacion inicial razonable

	L_VentanaImagen venCaja(0, 0, im.lx, im.ly, "Caja (click=salir)");

	rec.fijaReconocedor(&sift, &tree, &hough);
	rec.debugActivo(true);
	hough.hash.dxBin = 0.05;
	hough.hash.dyBin = 0.05;
	hough.hash.doctBin = 0.5;

	while(true)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		lins.resize(0);

		//while (venCaja.mouseE() != L_MousePresiona)
		//	Fl::check();
		if (venCaja.mouseE() != L_MouseLeido)
			break;

		venCaja.mouseE() = L_MouseLeido;

		imTmp = im;

		//desArr.generarDescriptorListaOriginal(sift.desL);
		rec.corrPixBN.verif = false;
		hough.fusion.activo = false;
		rec.debugActivo(0);
		//
		// Procesamiento para encontrar la mejor pose
		//desArr.generarDescriptorLista(sift.desL, pose);
		desArr.generarDescriptorLista_normalImprecisa(sift.desL, pose);
		rec.agregaImagenesEntren_usa_desL(sift.desL);
		rec.finalizaEntren();
		rec.procesaImagenesPrueba(imx);

		if (rec.traL.size() > 0)
		{
			arrRay.resize(rec.traL.root->c->calL.size());
			arrPos.resize(rec.traL.root->c->calL.size());

			i = 0;
			for (cal = rec.traL.root->c->calL.root; cal!=NULL; cal=cal->sig)
			{
				cam.calcRayo(arrRay[i], cal->c.dPru->x0, cal->c.dPru->y0);
				arrPos[i] = cal->c.dRef->d3d->cen;
				lins.drawLine((int)cal->c.dPru->x0, (int)cal->c.dPru->y0, (int)cal->c.dRef->x0, (int)cal->c.dRef->y0, 255, 255, 255);
				i++;
			}
			opt.lenVector = 7; // puntocuaternion
			opt.nIterationsMax = 8;  // ojala que converja con esto

			pose.calcPose_movido_a_fijo_Nrayos_sqrt(arrPos, arrRay, opt);

			p3p[0].fijaXYZ(0, 0, 0);
			p3p[1].fijaXYZ(0, -147*asp, 0);
			p3p[2].fijaXYZ(0, -147*asp, -167);
			p3p[3].fijaXYZ(0, -0, -167);

			p3p[0] = pose.movido_a_fijo(p3p[0]);
			p3p[1] = pose.movido_a_fijo(p3p[1]);
			p3p[2] = pose.movido_a_fijo(p3p[2]);
			p3p[3] = pose.movido_a_fijo(p3p[3]);

			for (i=0; i<4; i++)
			{
				ray.define(p3p[i].y / p3p[i].x, p3p[i].z / p3p[i].x);
				cam.calcPixel(ray, x0[i], y0[i]);
			}
			for (i=0; i<4; i++)
				lins.drawLine(x0[i], y0[i], x0[(i+1)%4], y0[(i+1)%4], 255, 255, 255);

			ptr = pose.calcPose3D();

			printf("x: %.2f,  y: %.2f,  z: %.2f\n", pose.pos.x, pose.pos.y, pose.pos.z);
			printf("pan: %.0f  tilt: %.0f  roll: %.0f\n", ptr.ori.pan*180/M_PI, ptr.ori.tilt*180/M_PI, ptr.ori.roll*180/M_PI);
		}
		//rec.desLPru.dibujaFlechas(linsCamara);
		for (int j=0; j<imTmp.ly; j++)
			for (int i=0; i<imTmp.lx; i++)
				for (int k=0; k<3; k++)
					imTmp.pix(i,j,k) = (L_uchar)(imTmp.pix(i,j,k)*0.4);
		imTmp.genDrawing(lins);
		venCaja.dibuja(imTmp);
		rec.reseteaPrueba();
		rec.reseteaEntren();
		lins.resize(0);

		venCaja.check();
	}
	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)


#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_dibModelo(const L_Descriptor3DArreglo &des3DArr, int xIni)
{
	L_VentanaImagen venCaja(0, 0, 640, 480, "Modelo (click=salir)");
	L_VentanaImagen venPanTilt(0, 500, 200, 200, "Pan - tilt");
	L_VentanaImagen venXY(250, 500, 200, 200, "x,y");
	L_CamaraPinhole cam;
	L_Pose3D pose;
	L_HomogeneousMatrix H;
	L_Pose3D_cuat posec;
	L_DescriptorLista desL;
	L_ShapeArray lins;
	L_ImageRGBUchar im;
	L_DescriptorNodo *desn;
	L_CoordsCart3D prom;
	L_Descriptor3DArreglo des3DA;
	int i, u;

	des3DA = des3DArr;

	im.reallocate(640, 480);
	cam.fijaParametrosPhillips(640, 480);

	for (i=0; i<(int)des3DArr.size(); i++)
		des3DA[i].cam = cam;

	for (i=0; i<(int)des3DArr.size(); i++)
		des3DA[i].sigma1 = sqrt(des3DA[i].sigma1)*10; // Para eitar que algunos valores se escapen

	pose.fijaCero();
	pose.pos.x = xIni;

	prom.fijaCero();
	for (i=0; i<(int)des3DA.size(); i++)
		prom += des3DA[i].cen;
	prom /= des3DA.size();

	printf("Centro del object: %g %g %g\n", prom.x, prom.y, prom.z);

	while(true)
	{
		if (venCaja.mouseE() != L_MouseLeido)
			break;
		if (venPanTilt.mouseE() != L_MouseLeido)
		{
			pose.ori.pan = (venPanTilt.mouseX() - 100)/100.0 * 2*M_PI;
			pose.ori.tilt = (100 - venPanTilt.mouseY())/100.0 * 2*M_PI;
			pose.imprime_gr("pose");
			venPanTilt.mouseE() = L_MouseLeido;
		}
		if (venXY.mouseE() != L_MouseLeido)
		{
			pose.pos.x = (100-venXY.mouseY()) * 5;
			pose.pos.y = (100-venXY.mouseX()) * 5;
			pose.imprime_gr("pose");
			venXY.mouseE() = L_MouseLeido;
		}
		H.fijaPose3D(pose);
		posec = H.calcPose3D_cuat();
		des3DA.generarDescriptorLista(desL, posec);
		u=0;
		for (desn=desL.root; desn!=NULL; desn=desn->sig)
		{
			//lins.drawArrow((int)(desn->c.x0), (int)(desn->c.y0), 20, desn->c.ang, (L_uchar)((u*10)%255), (L_uchar)((u*25)%255), (L_uchar)((u*100)%255));
			lins.drawArrow((int)(desn->c.x0), (int)(desn->c.y0), desn->c.sigma0*10, desn->c.ang, (L_uchar)((u*10)%255), (L_uchar)((u*25)%255), (L_uchar)((u*100)%255));
			u++;
		}
		im.setZero();
		im.genDrawing(lins);
		venCaja.dibuja(im);
		desL.clear();
		lins.resize(0);
		venCaja.check();
	}
	return 0;

}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)


#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camCalcPosePuntInt_1(int argc, char *argv[])
{
	int numPan, numTilt, i, j, x0, y0;
	double dpan, dtilt;
	double pan[16], tilt[16];
	bool entrenado;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;
	L_CapturadorImagen capt;
	L_ShapeArray linsCamara;
	L_Array<L_ImageRGBUchar> imRefArr;
	L_Array<const L_ImageGrayDouble *> imarr;

	L_OpenSURF sift;
	L_KdTreeBBF tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;
	
	printf("Ingrese el numero de angulos horizontales: ");
	while (scanf("%d", &numPan) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese el numero de angulos verticales: ");
	while(scanf("%d", &numTilt) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese la diferencia de angulo horizontal (en grados): ");
	while(scanf("%lf", &dpan) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese la diferencia de angulo vertical (en grados): ");
	while(scanf("%lf", &dtilt) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venVistas(0, im.ly+50, im.lx, im.ly, "Vistas (click=salir)");

	dpan = dpan * M_PI / 180;
	dtilt = dtilt * M_PI / 180;
	for (i=0; i<numPan; i++)
		pan[i] = dpan*i - dpan*((numPan-1.0)/2);
	for (j=0; j<numTilt; j++)
		tilt[j] = dtilt*j - dtilt*((numTilt-1.0)/2);

	imRefArr.resize(numPan*numTilt);

	printf("Ingrese las vistas generando rectangulos con el mouse\n");

	i=0;
	j=0;

	rec.fijaReconocedor(&sift, &tree, &hough);
	while(i < numPan && j < numTilt)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);

		if (venCam.mouseE() == L_MousePresiona)
		{
			venCam.mouseE() = L_MouseLeido;
			x0 = venCam.mouseX();
			y0 = venCam.mouseY();
			entrenado = false; // Para subir el frame rate
		}
		if (venCam.mouseE() == L_MouseArrastra)
		{
			venCam.mouseE() = L_MouseLeido;
			entrenado = false; // Para subir el frame rate
			linsCamara.drawRectangle(x0, y0, venCam.mouseX(), venCam.mouseY());
		}
		if (venCam.mouseE() == L_MouseSuelta)
		{
			venCam.mouseE() = L_MouseLeido;
			if (venCam.mouseX()-x0 > 0 && venCam.mouseY()-y0 > 0) // Por si hay problemas de sincronismo
			{
				int xM=venCam.mouseX()-x0+1, yM=venCam.mouseY()-y0+1;
				if (xM > capt.im.lx - x0)
					xM = capt.im.lx - x0 -1;
				if (yM > capt.im.ly - y0)
					yM = capt.im.ly - y0 -1;
				imTrozo.subImageOf(capt.im, x0, y0, xM, yM);
				imTrozox = imTrozo;
				rec.agregaImagenEntrenCopiaObjetoDe(imTrozox);
				entrenado = true;
				imRefArr[i+j*numPan] = imTrozo;
				printf("%d de %d\n", i+j*numPan+1, numPan*numTilt);
				i++;
				if (i >= numPan)
				{
					i = 0;
					j++;
				}
			}
		}
		if (venVistas.mouseE() != L_MouseLeido)
			return 0;
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibuja(imTmp);
		venCam.check();
	}
	rec.finalizaEntren();

	while(true)
	{
		if (venVistas.mouseE() != L_MouseLeido)
			return 0;

		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		venCam.dibuja(im);

		rec.procesaImagenesPrueba(imx);
		rec.dibujaTransformaciones(imTmp, 1, 100, 1, 1);

		venVistas.resize(imTmp.lx, imTmp.ly);
		venVistas.dibuja(imTmp);
		venCam.check();
	}

	printf("presione ENTER para continuar\n");
	getchar();

	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)


#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
int L_ProgramasRecPuntInt::main_camCalcPosePuntInt_2(int argc, char *argv[])
{
	int numPan, numTilt, i, j, x0, y0;
	double dpan, dtilt;
	double pan[16], tilt[16];
	bool entrenado;
	L_ImageRGBUchar im, imTrozo, imTmp;
	L_ImageGrayDouble imx, imTrozox;
	L_CapturadorImagen capt;
	L_ShapeArray lins, linsCamara;
	L_Array<L_ImageRGBUchar> imRefArr;
	L_Array<const L_ImageGrayDouble *> imarr;
	L_DescriptorLista desL;

	L_OpenSURF sift;
	L_Array<L_KdTreeBBF> tree;
	L_GenTransfAfinHough hough;
	L_Array<L_RecPuntInt> rec;
	
	printf("Ingrese el numero de angulos horizontales: ");
	while (scanf("%d", &numPan) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese el numero de angulos verticales: ");
	while (scanf("%d", &numTilt) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese la diferencia de angulo horizontal (en grados): ");
	while (scanf("%lf", &dpan) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	printf("Ingrese la diferencia de angulo vertical (en grados): ");
	while(scanf("%lf", &dtilt) != 1)
		{fflush(stdin); printf("Ingrese un numero seguido de ENTER:\n");}

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara disponible. Presione ENTER");
		getchar();
		return 1;
	}
	im = capt.im;

	L_VentanaImagen venCam(0, 0, im.lx, im.ly, "Hacer cuadrado con el mouse");
	L_VentanaImagen venVistas(0, im.ly+50, im.lx, im.ly, "Vistas (click=salir)");

	dpan = dpan * M_PI / 180;
	dtilt = dtilt * M_PI / 180;

	for (i=0; i<numPan; i++)
		pan[i] = dpan*i - dpan*((numPan-1.0)/2);
	for (j=0; j<numTilt; j++)
		tilt[j] = dtilt*j - dtilt*((numTilt-1.0)/2);

	imRefArr.resize(numPan*numTilt);

	rec.resize_swapping(numPan*numTilt);
	tree.resize_swapping(numPan*numTilt);
	imarr.resize(numPan*numTilt+1);

	printf("Ingrese las vistas generando rectangulos con el mouse\n");

	for (i=0; i<numPan; i++)
		for (j=0; j<numTilt; j++)
			rec[i+j*numPan].fijaReconocedor(&sift, &tree[i+j*numPan], &hough);

	i=0;
	j=0;

	while(i < numPan && j < numTilt)
	{
		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		linsCamara.resize(0);

		if (venCam.mouseE() == L_MousePresiona)
		{
			venCam.mouseE() = L_MouseLeido;
			x0 = venCam.mouseX();
			y0 = venCam.mouseY();
			entrenado = false; // Para subir el frame rate
		}
		if (venCam.mouseE() == L_MouseArrastra)
		{
			venCam.mouseE() = L_MouseLeido;
			entrenado = false; // Para subir el frame rate
			linsCamara.drawRectangle(x0, y0, venCam.mouseX(), venCam.mouseY());
		}
		if (venCam.mouseE() == L_MouseSuelta)
		{
			venCam.mouseE() = L_MouseLeido;
			if (venCam.mouseX()-x0 > 0 && venCam.mouseY()-y0 > 0) // Por si hay problemas de sincronismo
			{
				int xM=venCam.mouseX()-x0+1, yM=venCam.mouseY()-y0+1;
				if (xM > capt.im.lx - x0)
					xM = capt.im.lx - x0 -1;
				if (yM > capt.im.ly - y0)
					yM = capt.im.ly - y0 -1;
				imTrozo.subImageOf(capt.im, x0, y0, xM, yM);
				imTrozox = imTrozo;
				rec[i+j*numPan].agregaImagenEntrenCopiaObjetoDe(imTrozox);
				rec[i+j*numPan].finalizaEntren();
				entrenado = true;
				imRefArr[i+j*numPan] = imTrozo;
				printf("%d de %d\n", i+j*numPan+1, numPan*numTilt);
				i++;
				if (i >= numPan)
				{
					i = 0;
					j++;
				}
			}
		}
		if (venVistas.mouseE() != L_MouseLeido)
			return 0;
		imTmp = im;
		imTmp.genDrawing(linsCamara);
		venCam.dibuja(imTmp);
		venCam.check();
	}

	for (i=0; i<numPan; i++)
			for (j=0; j<numTilt; j++)
				imarr[i+j*numPan+1] = rec[i+j*numPan].imsref[0];

	while(true)
	{
		if (venVistas.mouseE() != L_MouseLeido)
			return 0;

		if (capt.capturarImagen() == false)
		{
			printf("No hay camara disponible. Presione ENTER");
			getchar();
			return 2;
		}
		im = capt.im;
		imx = im;

		venCam.dibuja(im);

		rec[0].reseteaPrueba();
		rec[0].procesaImagenesPrueba(imx);
		rec[0].desLPru.moveListTo(desL);

		for (i=0; i<numPan; i++)
		{
			for (j=0; j<numTilt; j++)
			{
				rec[i+j*numPan].reseteaPrueba();
				rec[i+j*numPan].procesaImagenesPrueba_usa_desL(desL, &imx);
				if (rec[i+j*numPan].traL.size() > 0)
				{
					rec[i+j*numPan].traL.root->c->calL.genLineas(lins);
					lins.setImageNumber_ref(i+j*numPan+1);
					lins.swap(linsCamara);
				}
				//rec[i+j*numPan].dibujaTransformaciones(imTmp, 1, 100, 1, 1);
				rec[i+j*numPan].desLPru.moveListTo(desL);
			}
		}
		imarr[0] = &imx;
		imTmp = im;
		imTmp.genDrawingMatches(imarr, linsCamara);
		venVistas.resize(imTmp.lx, imTmp.ly);
		venVistas.dibuja(imTmp);
		linsCamara.resize(0);
		desL.clear();
		venCam.check();
	}

	printf("presione ENTER para continuar\n");
	getchar();

	return 0;
}
#endif // defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)

