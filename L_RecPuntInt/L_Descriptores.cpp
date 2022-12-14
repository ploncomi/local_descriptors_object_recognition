#include "L_Descriptores.h"
#include "L_GenPuntInt.h"
//#include "../L_LoweSIFT/L_Lowe.h"
//#include "../OpenSURF/L_OpenSURF.h"

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


L_SIFTGenFnDescriptorArreglo L_SIFTGenFnArr(true); // El unico object de tipo L_SIFTGenFnDescriptorArreglo

bool L_DescriptorLista::leeArchivoPato_minu(const char *nomarch, L_ImageGrayDouble &imDummy, L_PuntIntLista &punL, bool agregaEspejo)
{
	FILE *fp;
	int ene, v;
	int tipo;
	char str[10];
	double ang;
	fp=fopen(nomarch,"r");
	if (fp==NULL)
		return false;
	if (fgets(str, 10, fp) == NULL)
		{fclose(fp); printf("L_DescriptorLista::leeArchivoPato_minu() : formato de archivo incorrecto");return false;}
	if (strcmp(str,"MINUPATO\n")!=0)
	{
		fclose(fp);
		return false;
	}
	if (fscanf(fp, "%d%d", &v, &ene) != 2)
		{fclose(fp); printf("L_DescriptorLista::leeArchivoPato_minu() : formato de archivo incorrecto");return false;}
	if (v!=0)
	{
		printf("Error en formato de archivo en L_DireccLista::leeArchivoPato()\n");
		return false;
	}
	for (; ene>0; ene--)
	{
		*pult=new L_DireccNodo;
		*punL.pult=new L_PuntIntNodo;
		if (fscanf(fp, "%lf%lf%lf%lf%d", &(*punL.pult)->c.x0, &(*punL.pult)->c.y0, &(*punL.pult)->c.sigma0, &ang, &tipo) != 5)
			{fclose(fp); printf("L_DescriptorLista::leeArchivoPato_minu() : formato de archivo incorrecto");return false;}
		(*punL.pult)->c.tipoP=MINUCIA;
		(*punL.pult)->c.intens = 1.0;
		(*punL.pult)->c.signo = (tipo==0) ? 1 : -1; // Termino > 0,  bifurcacion < 0
		(*punL.pult)->c.imorig=&imDummy;
		(*punL.pult)->c.imorig_lx=imDummy.lx;
		(*punL.pult)->c.imorig_ly=imDummy.ly;
		(*pult)=*punL.pult;
		(*pult)->c.ang=ang;
		if (agregaEspejo)
		{
			ang =(*pult)->c.ang+M_PI;
			pult=&(*pult)->sig;
			n++;
			*pult=new L_DireccNodo;
			(*pult)=*punL.pult;
			(*pult)->c.ang=ang;
		}
		pult=&(*pult)->sig;
		n++;
		punL.pult=&(*punL.pult)->sig;
		punL.n++;
	}
	fclose(fp);
	return true;
}

bool L_DescriptorLista::grabaArchivoPato_minu(const char *nomarch)
{
	FILE *fp;
	L_DireccNodo *des;
	fp=fopen(nomarch,"w");
	if (fp==NULL)
		return false;
	fputs("MINUPATO\n",fp);
	if (root==NULL)
	{
		fprintf(fp, "0\n0\n");
		fclose(fp);
		return true;
	}
	fprintf(fp,"%d\n%d\n", 0, n);
	for(des=root; des!=NULL; des=des->sig)
	{
		fprintf(fp, "%lf %lf %lf %lf %d", des->c.x0, des->c.y0, des->c.sigma0, des->c.ang, (des->c.intens > 0));
		fprintf(fp,"\n");
	}
	fclose(fp);
	return true;
}



void L_Descriptor::creaVectorPromedioDe(const L_Descriptor &d1, const L_Descriptor &d2)
{
	int i;
	if (d1.nt != d2.nt || d1.vector == NULL || d2.vector == NULL)
	{
		printf("Error en argumentos de L_Descriptor::creaVectorPromedioDe()\n");
		return;
	}
	*this = d1;
	for (i=0; i<d1.nt; i++)
		vector[i]=(d1.vector[i] + d2.vector[i])*0.5;
}

double L_Descriptor::distEuclidCuadr(const L_Descriptor & other) const
{
	int i;
	double sumacuad=0;
	throw_L_ArgException_if(nt!=other.nt, "L_Descriptor::distEuclidCuadr");
	for (i=0; i<nt; i++)
		sumacuad+=(vector[i]-other.vector[i])*(vector[i]-other.vector[i]);
	return sumacuad;
}

double L_Descriptor::distEuclidCuadr_mismoTipo(const L_Descriptor & other) const
{
	int i;
	double sumacuad=0;
	throw_L_ArgException_if(nt!=other.nt, "L_Descriptor::distEuclidCuadr");
	for (i=0; i<nt; i++)
		sumacuad+=(vector[i]-other.vector[i])*(vector[i]-other.vector[i]);
	if (tipoP != other.tipoP || tipoD != other.tipoD)
		sumacuad += 1000;
	return sumacuad;
}

bool L_Descriptor::normalizaHist()
{
	double *v1, *v1f;
	double sum=0;
	v1=vector;
	v1f=vector+nt-1;
	while (v1<=v1f)
		sum+=*(v1++);
	if (sum==0)
		return false;
	v1=vector;
	while (v1<=v1f)
		*(v1++)/=sum;
	return true;
}

bool L_Descriptor::normalizaHistCota(double maxValorComponente)
{
	double *v1, *v1f;
	if (normalizaHist()==false)
		return false;
	v1=vector;
	v1f=vector+nt-1;
	while (v1<=v1f)
		if (*v1>maxValorComponente)
			*v1=maxValorComponente;
	return normalizaHist();
}

int L_Descriptor::cmpPosSigmaAng(const L_Descriptor &other) const
{
#define L_D_CPSA_CMP(X) {if ((X) > other.X) return 1; else if ((X) < other.X) return -1;}
	L_D_CPSA_CMP(x0);
	L_D_CPSA_CMP(y0);
	L_D_CPSA_CMP(sigma0);
	L_D_CPSA_CMP(ang);
	return 0;
#undef L_D_CPSA_CMP
}

int L_Descriptor::cmpPtr(const L_Descriptor &other) const
{
	if (this > &other)
		return 1;
	else if (this < &other)
		return -1;
	return 0;
}

int L_Descriptor::cmpMula(const L_Descriptor &other) const
{
	int i;
	for (i=0; i<nt; i++)
	{
		if (vector[i]>other.vector[i])
			return 1;
		else if (vector[i]<other.vector[i])
			return -1;
	}
	return 0;
}

int L_Descriptor::qCmpMula_nodo(const void *pp1, const void *pp2)
{
	return ((**(L_DescriptorNodo **)pp1).c).cmpMula((**(L_DescriptorNodo **)pp2).c);
}

int L_Descriptor::qCmpMulaPtr_nodo(const void *pp1, const void *pp2)
{
	return ((**(L_DescriptorPtrNodo**)pp1)).c.ptr->c.cmpMula((**(L_DescriptorPtrNodo**)pp2).c.ptr->c);
}

int L_Descriptor::qCmpIntens(const void *p1, const void *p2)
{
	return
		( (*(L_Descriptor*)p1).intens > (*(L_Descriptor*)p2).intens )
		- ( (*(L_Descriptor*)p1).intens < (*(L_Descriptor*)p2).intens );
}

int L_Descriptor::qCmpIntens_nodo(const void *pp1, const void *pp2)
{
	return
		( (**(L_DescriptorNodo**)pp1).c.intens > (**(L_DescriptorNodo**)pp2).c.intens )
		- ( (**(L_DescriptorNodo**)pp1).c.intens < (**(L_DescriptorNodo**)pp2).c.intens );
}

int L_Descriptor::qCmpIntensInv_nodo(const void *pp1, const void *pp2)
{
	return
		( (**(L_DescriptorNodo**)pp1).c.intens < (**(L_DescriptorNodo**)pp2).c.intens )
		- ( (**(L_DescriptorNodo**)pp1).c.intens > (**(L_DescriptorNodo**)pp2).c.intens );
}


double L_Descriptor::medidaHarris(L_Array<L_ImageGrayDouble> &imArrTmp, double sd, double si, double alfa)
{
	if (imorig != NULL)  // Pueden ser descriptores sin referencia a la imagen original
		return L_ImageGrayDouble::harrisIndividual(imArrTmp, *imorig, x0, y0, sd, si, alfa);
	else
		return 0.001; // Un cornerness grande para que no se descarte
}


void L_Descriptor::computeCovariance3DFromDisparity(double sigmaPixNormX, double sigmaPixNormY, double sigmaPixNormDisparity, double baselinemm)
{
	double u, v, d, sigma_depth, covu, covv, covd;

	u = ray.tanIzq;
	v = ray.tanArr;
	d = xAde;

	sigma_depth = sigmaPixNormDisparity / baselinemm * d * d;

	covu = sigmaPixNormX * sigmaPixNormX;
	covv = sigmaPixNormY * sigmaPixNormY;
	covd = sigma_depth * sigma_depth;

	cov[0][0] = covd;
	cov[0][1] = u * covd;
	cov[0][2] = v * covd;
	cov[1][0] = u * covd;
	cov[1][1] = d * d * covu + u * u * covd;
	cov[1][2] = u * covd * v;
	cov[2][0] = v * covd;
	cov[2][1] = u * covd * v;
	cov[2][2] = d * d * covv + v * v * covd;
}


void L_Descriptor::writeElementTxt(FILE *fp)
{
	int i;
	fprintf(fp, "%.2f %.2f %lf %d %d %lf %d %d\n", x0, y0, sigma0, tipoP, nt, ang, nobj, tipoD);
	for (i=0; i<nt; i++)
		fprintf(fp, "%lf ", vector[i]);
	fprintf(fp,"\n");
}

void L_Descriptor::readElementTxt(FILE *fp, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen)
{
	int i;
	int i_tipoP=(int)tipoP;
	int i_tipoD=(int)tipoD;
	if (fscanf(fp, "%lf%lf%lf%d%d%lf%d%d", &x0, &y0, &sigma0, &i_tipoP, &nt, &ang, &nobj, &i_tipoD) != 8)
		{printf("L_Descriptor::readElementTxt() : formato de archivo incorrecto\n"); return;}
	for (i=0; i<nt; i++)
		if (fscanf(fp, "%lf", &vector[i]) != 1)
			{printf("L_Descriptor::readElementTxt() : formato de archivo incorrecto\n"); return;}
	imorig=im;
	refGen=gen;
}

void L_Descriptor::writeElementBin(FILE *fp)
{
	fwrite(this, sizeof(*this), 1, fp);
	//fwrite(vector, sizeof(vector[0]), nt, fp);
}

void L_Descriptor::readElementBin(FILE *fp, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen)
{
	if (fread(this, sizeof(*this), 1, fp) != 1)
		{printf("L_Descriptor::readElementBin() : formato de archivo incorrecto\n"); return;}
	//fread(vector,sizeof(vector[0]), nt, fp);
	imorig=im;
	if (gen->tipoP == P_UNION && gen->tipoD == D_UNION) // gen es un arreglo de generadores, buscar cual corresponde
	{
		L_GenDescrArr *g = (L_GenDescrArr *)gen;
		int i;
		for (i=0; i<g->n; i++)
		{
			if (g->arr[i]->tipoP == tipoP && g->arr[i]->tipoD == tipoD)
			{
				refGen = g->arr[i];
				break;
			}
		}
		throw_L_ArgException_if(i == g->n, "L_Descriptor::readElementBin(.,L_GenDescrArr)");
	}
	else
		refGen=gen;
	nobj = refGen->nobj;
}

void L_Descriptor::readElementBin(FILE *fp, va_list va)
{
	L_ImageGrayDouble *im = va_arg(va, L_ImageGrayDouble *);
	L_GenDescrPuntInt *gpi = va_arg(va, L_GenDescrPuntInt *);
	readElementBin(fp, im, gpi);
}

int L_Descriptor::writeElementMem(std::vector<char> &buf)
{
	int i0 = (int)buf.size();
	buf.resize(buf.size()+sizeof(L_Descriptor));// + sizeof(double)*nt;
	memcpy(&(buf[0])+i0, this, sizeof(L_Descriptor)); // OJO
	return sizeof(L_Descriptor);// + sizeof(double)*nt;
}

int L_Descriptor::readElementMem(const std::vector<char> &buf, int i0, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen)
{
	memcpy(this, &buf[i0], sizeof(*this));
	//memcpy(vector, &buf[i0]+sizeof(L_Descriptor), sizeof(double)*nt);
	imorig=im;
	if (gen->tipoP == P_UNION && gen->tipoD == D_UNION) // gen es un arreglo de generadores, buscar cual corresponde
	{
		L_GenDescrArr *g = (L_GenDescrArr *)gen;
		int i;
		for (i=0; i<g->n; i++)
		{
			if (g->arr[i]->tipoP == tipoP && g->arr[i]->tipoD == tipoD)
			{
				refGen = g->arr[i];
				break;
			}
		}
		throw_L_ArgException_if(i == g->n, "L_Descriptor::readElementBin(.,L_GenDescrArr)");
	}
	else
		refGen=gen;
	return sizeof(L_Descriptor);// + sizeof(double)*nt;
}

int L_Descriptor::readElementMem(const std::vector<char> &buf, int i0, va_list va)
{
	const L_ImageGrayDouble *im = va_arg(va, L_ImageGrayDouble *);
	L_GenDescrPuntInt *gpi = va_arg(va, L_GenDescrPuntInt *);
	return readElementMem(buf, i0, im, gpi);
}

L_Descriptor L_Descriptor::creaDescriptorAzar(int lx, int ly, int dim, L_GenDescrPuntInt *gen)
{
	L_Descriptor d;
	double azar;
	int i;
	d.x0 = L_RANDOM(0,lx);
	d.y0 = L_RANDOM(0,ly);
	d.ang = L_RANDOM(0, 2*M_PI);
	azar = L_RANDOM(0, 3);
	d.sigma0 = pow(2.0,azar); // Asi se concentra en los numeros mas bajos
	d.nt = dim;
	for (i=0; i<dim; i++)
		d.vector[i] = L_RANDOM(-1, 1);
	d.intens = L_RANDOM(1.0e-3, 1.0e-6);
	d.hasDepth = 0;

	d.normalizaHist();
	d.imorig_lx = lx;
	d.imorig_ly = ly;
	d.imorig = NULL;
	d.cam = 0; // Por defecto
	d.pir_x = d.x0 / d.sigma0;
	d.pir_y = d.y0 / d.sigma0;
	d.pir_n = 1;
	d.refGen = gen;
	d.tipoP = gen->tipoP;
	d.tipoD = gen->tipoD;
	return d;
}

L_Descriptor L_Descriptor::creaDescriptorAzarSemilla(int lx, int ly, int dim, L_GenDescrPuntInt *gen, L_Rand &randF)
{
	L_Descriptor d;
	double azar;
	int i;
	d.x0 = randF.random(0,lx);
	d.y0 = randF.random(0,ly);
	d.ang = randF.random(0, 2*M_PI);
	azar = randF.random(0, 3);
	d.sigma0 = pow(2.0,azar); // Asi se concentra en los numeros mas bajos
	d.nt = dim;
	for (i=0; i<dim; i++)
		d.vector[i] = randF.random(-1, 1);
	d.intens = randF.random(1.0e-3, 1.0e-6);
	d.hasDepth = 0;

	d.normalizaHist();
	d.imorig_lx = lx;
	d.imorig_ly = ly;
	d.imorig = NULL;
	d.cam = 0; // Por defecto
	d.pir_x = d.x0 / d.sigma0;
	d.pir_y = d.y0 / d.sigma0;
	d.pir_n = 1;
	d.refGen = gen;
	d.tipoP = gen->tipoP;
	d.tipoD = gen->tipoD;
	return d;
}

#ifdef IMGFEATURES_H
void L_Descriptor::OP_assign(const feature &feat)
{
	int i;
	x0 = feat.x;
	y0 = feat.y;
	sigma0 = feat.scl;
	ang = feat.ori;
	nt = feat.d;
	refGen = NULL;
	imorig = NULL;
	for (i=0; i<nt; i++)
		vector[i] = feat.descr[i];
}
#endif // IMGFEATURES_H

#ifdef IMGFEATURES_H
void L_Descriptor::copyTo(feature &feat)
{
	int i;
	feat.x = x0;
	feat.y = y0;
	feat.scl = sigma0;
	feat.ori = ang;
	feat.d = nt;
	for (i=0; i<nt; i++)
		feat.descr[i] = vector[i];
}
#endif // IMGFEATURES_H

void L_DescriptorLista::dibujaFlechas(L_ShapeArray &lins, int r, int g, int b, int grosor)
{
	L_DescriptorNodo *ptr;
	double xF, yF;
	double xI, yI;
	double xD, yD;
	int n1, n2, n3;
	L_uchar R, G, B;
	const double fact=3.0;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		n1=(int)(ptr->c.ang*16/(2*M_PI))%16*16;
		n2=(int)( ptr->c.x0/5 );
		n3=(int)( ptr->c.y0/5 );
		if (r==-1)
			R=(L_uchar)(5*n1+2*n2+n3);
		else
			R=(L_uchar)r;
		if (g==-1)
			G=(L_uchar)(2*n1+n2+5*n3+128);
		else
			G=(L_uchar)g;
		if (b==-1)
			B=255-(R+G)/2;
		else
			B=(L_uchar)b;
		xF=+1.0*fact*ptr->c.sigma0*cos(ptr->c.ang);
		yF=-1.0*fact*ptr->c.sigma0*sin(ptr->c.ang);
		xI=+0.6*fact*ptr->c.sigma0*cos(ptr->c.ang+20.0*M_PI/180);
		yI=-0.6*fact*ptr->c.sigma0*sin(ptr->c.ang+20.0*M_PI/180);
		xD=+0.6*fact*ptr->c.sigma0*cos(ptr->c.ang-20.0*M_PI/180);
		yD=-0.6*fact*ptr->c.sigma0*sin(ptr->c.ang-20.0*M_PI/180);
		if (grosor==1)
		{
			// linea principal de la flecha
			lins.drawLine((int)(ptr->c.x0+0.5),(int)(ptr->c.y0+0.5),(int)(ptr->c.x0+0.5+xF),(int)(ptr->c.y0+0.5+yF), R,G,B);
			// otras dos lineas de la flecha
			lins.drawLine((int)(ptr->c.x0+0.5+xF),(int)(ptr->c.y0+0.5+yF), (int)(ptr->c.x0+0.5+xD), (int)(ptr->c.y0+0.5+yD), R,G,B);
			lins.drawLine((int)(ptr->c.x0+0.5+xF),(int)(ptr->c.y0+0.5+yF), (int)(ptr->c.x0+0.5+xI), (int)(ptr->c.y0+0.5+yI), R,G,B);
		}
		else
		{
			int i, j;
			for (i=0; i<grosor; i++)
			{
				for (j=1; j<grosor; j++)
				{
					// linea principal de la flecha
					lins.drawLine((int)(ptr->c.x0+0.5)+i,(int)(ptr->c.y0+0.5)+j,(int)(ptr->c.x0+0.5+xF)+i,(int)(ptr->c.y0+0.5+yF)+j, R,G,B);
					// otras dos lineas de la flecha
					lins.drawLine((int)(ptr->c.x0+0.5+xF)+i,(int)(ptr->c.y0+0.5+yF)+j, (int)(ptr->c.x0+0.5+xD)+i, (int)(ptr->c.y0+0.5+yD)+j, R,G,B);
					lins.drawLine((int)(ptr->c.x0+0.5+xF)+i,(int)(ptr->c.y0+0.5+yF)+j, (int)(ptr->c.x0+0.5+xI)+i, (int)(ptr->c.y0+0.5+yI)+j, R,G,B);
				}
			}
		}
	}
}

void L_DescriptorLista::dibujaCuadrados(L_ShapeArray &lins)
{
	L_DescriptorNodo *ptr;
	double xD, yD;
	L_uchar R, G, B;
	const double fact=1.0;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		L_ShapeArray::genRandomColor(R, G, B);
		xD=+fact*ptr->c.sigma0*cos(ptr->c.ang);
		yD=-fact*ptr->c.sigma0*sin(ptr->c.ang);
		lins.drawLine((int)(ptr->c.x0+0.5+xD),(int)(ptr->c.y0+0.5+yD),(int)(ptr->c.x0+0.5+yD),(int)(ptr->c.y0+0.5-xD), R,G,B);
		lins.drawLine((int)(ptr->c.x0+0.5+yD),(int)(ptr->c.y0+0.5-xD),(int)(ptr->c.x0+0.5-xD),(int)(ptr->c.y0+0.5-yD), R,G,B);
		lins.drawLine((int)(ptr->c.x0+0.5-xD),(int)(ptr->c.y0+0.5-yD),(int)(ptr->c.x0+0.5-yD),(int)(ptr->c.y0+0.5+xD), R,G,B);
		lins.drawLine((int)(ptr->c.x0+0.5-yD),(int)(ptr->c.y0+0.5+xD),(int)(ptr->c.x0+0.5+xD),(int)(ptr->c.y0+0.5+yD), R,G,B);
	}
}

bool L_DescriptorLista::leeArchivoPato(const char *nomarch, const L_ImageGrayDouble &im, L_GenDescrPuntInt *gen, int nobj)
{
	FILE *fp;
	int i, n, v;
	char str[10];
	fp=fopen(nomarch,"r");
	if (fp==NULL)
		return false;
	if (fgets(str, 9, fp) == NULL)
		{printf("L_DescriptorLista::leeArchivoPato() : formato de archivo incorrecto\n"); fclose(fp); return false;}
	if (strcmp(str,"DESCPATO")!=0)
	{
		fclose(fp);
		return false;
	}
	if (fscanf(fp, "%d%d", &v, &n) != 2)
		{printf("L_DescriptorLista::leeArchivoPato() : formato de archivo incorrecto\n"); fclose(fp); return false;}
	for (; n>0; n--)
	{
		*pult=new L_DescriptorNodo;
		if (fscanf(fp, "%lf%lf%lf%lf", &(*pult)->c.x0, &(*pult)->c.y0, &(*pult)->c.sigma0, &(*pult)->c.ang) != 4)
			{printf("L_DescriptorLista::leeArchivoPato() : formato de archivo incorrecto\n"); fclose(fp); return false;}
		(*pult)->c.nobj=nobj;
		(*pult)->c.nt=v;
		(*pult)->c.imorig=&im;
		(*pult)->c.imorig_lx=im.lx;
		(*pult)->c.imorig_ly=im.ly;
		(*pult)->c.refGen=gen;
		(*pult)->c.tipoP = gen->tipoP;
		(*pult)->c.tipoD = gen->tipoD;
		(*pult)->c.intens = 1;
		for (i=0; i<v; i++)
			if (fscanf(fp, "%lf", &(*pult)->c.vector[i]) != 1)
				{printf("L_DescriptorLista::leeArchivoPato() : formato de archivo incorrecto\n"); fclose(fp); return false;}
		pult=&(*pult)->sig;
		this->n++;
	}
	*pult=NULL;
	fclose(fp);
	return true;
}

bool L_DescriptorLista::grabaArchivoPato(const char *nomarch)
{
	FILE *fp;
	int i;
	L_DescriptorNodo *des;
	fp=fopen(nomarch,"w");
	if (fp==NULL)
		return false;
	fputs("DESCPATO\n",fp);
	if (root==NULL)
	{
		fprintf(fp, "0\n0\n");
		fclose(fp);
		return true;
	}
	fprintf(fp,"%d\n%d\n", root->c.nt, n);
	for(des=root; des!=NULL; des=des->sig)
	{
		fprintf(fp, "%lf %lf %lf %lf", des->c.x0, des->c.y0, des->c.sigma0, des->c.ang);
		for (i=0; i<des->c.nt; i++)
			fprintf(fp," %lf", des->c.vector[i]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return true;
}

void L_DescriptorLista::transformaEnMatriz(L_Matrix &uv)
{
	int i=0;
	L_DescriptorNodo *p;
	uv.reallocate(n, 2);
	for (p=root; p!=NULL; p=p->sig, i++)
	{
		uv(i,0)=p->c.x0;
		uv(i,1)=p->c.y0;
	}
}

bool L_DescriptorLista::revisaCoherencia()
{
	double sx=0, sy=0, sv=0;
	int ene=0;
	L_DescriptorNodo *d;
	for(d=root; d!=NULL; d=d->sig)
	{
		sx += d->c.x0 * d->c.x0;
		sy += d->c.y0 * d->c.y0;
		sv += d->c.vector[0] * d->c.vector[0];
		ene++;
	}
	throw_L_ArgException_if(ene != n, "L_DescriptorLista::revisaCoherencia() : n incorrecto");
	return sx==sx && sy==sy && sv==sv && (sx+sy+sv) != 0; // Ver que no haya NaNs y que no contenga sólo ceros
}


void L_DescriptorLista::copia_a_arreglo(std::vector<L_Descriptor> &arr)
{
	int i;
	L_DescriptorNodo *p;
	i = (int)arr.size();
	arr.resize(arr.size() + size());
	for (p=root; p!=NULL; p=p->sig)
		arr[i++] = p->c;
}

void L_DescriptorLista::copia_a_arreglo_punteros(std::vector<L_Descriptor *> &arr)
{
	int i;
	L_DescriptorNodo *p;
	i = (int)arr.size();
	arr.resize(arr.size() + size());
	for (p=root; p!=NULL; p=p->sig)
		arr[i++] = &p->c;
}

void L_DescriptorLista::agregarRuidoUniforme(double deltaPix)
{
	L_DescriptorNodo *p;
	for (p=root; p!=NULL; p=p->sig)
	{
		p->c.x0 += L_RANDOM(-deltaPix,deltaPix);
		p->c.y0 += L_RANDOM(-deltaPix,deltaPix);
	}
}

void L_DescriptorLista::agregarRuidoUniformeSemilla(double deltaPix, L_Rand &randF)
{
	L_DescriptorNodo *p;
	for (p=root; p!=NULL; p=p->sig)
	{
		p->c.x0 += randF.random(-deltaPix,deltaPix);
		p->c.y0 += randF.random(-deltaPix,deltaPix);
	}
}

void L_DescriptorLista::eliminaLineasHarris(double cornernessMin, double sd, double si, double alfa)
{
	int r = L_GaussianHalfWidth(si);
	int an = 2*r+1;
	int lx, ly, nop_indiv, nop_ant; // Numero de multiplicaciones a hacer con cada opcion

	if (root == NULL || n == 0)
		return;

	lx = root->c.imorig_lx;
	ly = root->c.imorig_ly;

	//        convoluciones   gradientes y productos punto a punto
	nop_ant = lx*ly*2*an*4 + lx*ly*5;
	//           convol "sd"   convol "si"  gradientes y productos punto a punto
	nop_indiv = 2*an*an*an*n + 2*an*an*n + an*an*5*n;

	if (nop_ant < nop_indiv)
		eliminaLineasHarris_ant(cornernessMin, sd, si, alfa);
	else
		eliminaLineasHarris_indiv(cornernessMin, sd, si, alfa);
}

void L_DescriptorLista::eliminaLineasHarris_indiv(double cornernessMin, double sd, double si, double alfa)
{
	if (n <= 0)
		return;
	std::vector<double> x0(n);
	std::vector<double> y0(n);
	L_Array<bool> eliminar(n);
	int i;
	L_DescriptorNodo *des;
	const L_ImageGrayDouble *act = NULL;
	L_ImageGrayDouble imHarris;

	for(i=0, des=root; i<n && des!=NULL; i++, des=des->sig)
	{
		eliminar[i] = false;
		if (des->c.imorig==NULL || des->c.imorig->data() == NULL)
			continue;
		throw_L_ArgException_if(des->c.imorig != act && act!=NULL, "L_DescriptorLista::eliminaLineasHarris_indiv() : descriptores con distintas imagenes origen");		
		throw_L_ArgException_if(des->c.imorig == NULL, "L_DescriptorLista::eliminaLineasHarris_indiv() : imorig no definido");		
		act = des->c.imorig;
		if (des->c.x0 < 0 || des->c.y0 < 0 || des->c.x0 > act->lx || des->c.y0 > act->ly)
			continue;

		x0[i] = des->c.x0;
		y0[i] = des->c.y0;
	}
	imHarris.harrisEliminaLineas(cornernessMin, *act, &(x0[0]), &(y0[0]), size(), eliminar.data());
	erase_preserving_order(eliminar);
}


void L_DescriptorLista::eliminaLineasHarris_ant(double cornernessMin, double sd, double si, double alfa)
{
	if (n <= 0)
		return;
	L_Array<bool> eliminar(n);
	int i;
	int x0, y0;
	L_DescriptorNodo *des;
	const L_ImageGrayDouble *act = NULL;
	L_ImageGrayDouble imHarris;

	for(i=0, des=root; i<n && des!=NULL; i++, des=des->sig)
	{
		eliminar[i] = false;
		if (des->c.imorig==NULL || des->c.imorig->data() == NULL)
			continue;
		if (des->c.imorig != act)
		{
			act = des->c.imorig;
			imHarris.harris(*act, sd, si, alfa);
		}
		x0 = (int)des->c.x0;
		y0 = (int)des->c.y0;
		if (x0 < 0)
			x0 = 0;
		if (x0 >= act->lx)
			x0 = act->lx-1;
		if (y0 < 0)
			y0 = 0;
		if (y0 >= act->ly)
			y0 = act->ly-1;
		if (imHarris.pix(x0,y0) < cornernessMin)
			eliminar[i] = true;
	}
	erase_preserving_order(eliminar);
}


int L_DescriptorLista::cuentaVisibles(int lx, int ly, double borde)
{
	if (n <= 0)
		return 0;
	L_DescriptorNodo *des;
	int visibles = 0;
	for(des=root; des!=NULL; des=des->sig)
	{
		if (des->c.x0 > borde && des->c.y0 > borde && des->c.x0 <= lx-borde && des->c.y0 <= ly-borde && des->c.xAde > 1.0e-30)
			visibles++;
	}
	return visibles;
}


void L_DescriptorLista::eliminaNoVisibles(int lx, int ly, double borde)
{
	if (n <= 0)
		return;
	L_Array<bool> eliminar(n);
	int i;
	L_DescriptorNodo *des;

	for (i=0; i<n; i++)
		eliminar[i] = false;

	for(i=0, des=root; i<n && des!=NULL; i++, des=des->sig)
	{
		if (des->c.x0 < borde || des->c.y0 < borde || des->c.x0 >= lx-borde || des->c.y0 >= ly-borde || des->c.xAde < 1.0e-30)
			eliminar[i] = true;
	}
	erase_preserving_order(eliminar);
}



void L_DescriptorLista::pedirMatrizDatos(L_Matrix &xy)
{
	int i;
	L_DescriptorNodo *des;
	xy.reallocate(n, 2);
	for (i=0,des=root; i<n; i++,des=des->sig)
	{
		xy(i,0) = des->c.x0;
		xy(i,1) = des->c.y0;
	}
}

void L_DescriptorLista::pedirMatrizDatosRayos(L_Matrix &xy, const L_CamaraPinhole &cam)
{
	int i;
	L_DescriptorNodo *des;
	L_RayoCamara r;
	xy.reallocate(n, 2);
	for (i=0,des=root; i<n; i++,des=des->sig)
	{
		cam.calcRayo(r, des->c.x0, des->c.y0);
		xy(i,0) = r.tanIzq;
		xy(i,1) = r.tanArr;
	}
}


int L_DescriptorLista::eliminaPose_Ransac3Puntos(const L_CamaraPinhole &cam, double errMax, int nIter, int nMin, int nAceptAutom, const std::vector<L_CoordsCart3D> &pun, const std::vector<double> &puntajes, L_Pose3D_cuat &pose)
{
	L_TransfPoseProy2D estimadorPose;
	L_Matrix xy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	int res;
	int i;

	pedirMatrizDatosRayos(xy, cam);

	estimadorPose.ray.resize(xy.li);

	for (i=0; i<xy.li; i++)
	{
		estimadorPose.ray[i].tanIzq = xy(i,0);
		estimadorPose.ray[i].tanArr = xy(i,1);
	}

	estimadorPose.pun = pun;

	res = estimadorPose.calcProsacSimple(puntajes, seleccionados, errMax, nIter, nMin, nAceptAutom);
	if (res < nMin)
		return nMin;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	pose = estimadorPose.pose;
	return res;
}

void L_DescriptorLista::fillDepth(L_ImageGrayDouble &imDepth)
{
	L_DescriptorNodo *d;
	int u, v;
	for (d=root; d!=NULL; d=d->sig)
	{
		d->c.hasDepth = 0;
		u = (int)d->c.x0;
		v = (int)d->c.y0;
		if (u>0 && u<imDepth.cols()-1 && v>0 && v<imDepth.rows()-1 && imDepth.pix(u,v) > 0)
		{
			d->c.xAde = imDepth.pix(u,v);
			d->c.hasDepth = 1;
		}
	}
}

void L_DescriptorLista::fillCovarianceInverseDistance(double sigmaPixNormX, double sigmaPixNormY, double sigmaPixNormDisparity, double baselinemm)
{
	L_DescriptorNodo *n;

	for (n=root; n!=NULL; n=n->sig)
	{
		if (n->c.hasDepth == 0)
			continue;

		n->c.computeCovariance3DFromDisparity(sigmaPixNormX, sigmaPixNormY, sigmaPixNormDisparity, baselinemm);
	}
}

void L_DescriptorLista::fillCovarianceInverseDistance_v0(double sigmar0, double sigmaw)
{
	L_DescriptorNodo *n;
	double C[3][3], R[3][3], RD[3][3];
	double u, v, d, r, ru;
	for (n=root; n!=NULL; n=n->sig)
	{
		u = n->c.ray.tanIzq;
		v = n->c.ray.tanArr;
		d = sqrt(1 + u*u + v*v);
		ru = 1 / sqrt(1+u*u);

		R[0][0] = 1/d;
		R[1][0] = u/d;
		R[2][0] = v/d;

		R[0][1] = -u*ru;
		R[1][1] = ru;
		R[2][1] = 0;

		// x2 = y0*z1 - z0*y1;
		// y2 = z0*x1 - x0*z1;
		// z2 = x0*y1 - y0*x1;

		R[0][2] = R[1][0]*R[2][1] - R[2][0]*R[1][1];
		R[1][2] = R[2][0]*R[0][1] - R[0][0]*R[2][1];
		R[2][2] = R[0][0]*R[1][1] - R[1][0]*R[0][1];

		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				C[i][j] = 0;
		
		r = n->c.xAde;

		C[0][0] = sigmar0*sigmar0/(r*r);
		C[1][1] = sigmaw*sigmaw;
		C[2][2] = sigmaw*sigmaw;


		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				RD[i][j] = 0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				for (int k=0; j<3; k++)
					RD[i][j] += R[i][k]*C[k][j];  // RD = R*C

		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				n->c.cov[i][j] = 0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				for (int k=0; j<3; k++)
					n->c.cov[i][j] += RD[i][k]*C[j][k];  // cov = RD * CT
	}
}

#ifdef IMGFEATURES_H
void L_DescriptorLista::OP_assign(const feature *feat, int n)
{
	L_Descriptor d;
	int i;
	for (i=0; i<n; i++)
	{
		d.OP_assign(feat[i]);
		push_back(d);
	}
}
#endif // IMGFEATURES_H

#ifdef IMGFEATURES_H
void L_DescriptorLista::copyTo(feature *feat)
{
	L_DescriptorNodo *d;
	int i=0;
	for (d=root; d!=NULL; d=d->sig)
		d->c.copyTo(feat[i++]);
}
#endif // IMGFEATURES_H

int L_DescriptorPtr::qCmpMula(const void *p1, const void *p2)
{
	return (*(L_DescriptorPtrNodo **)p1)->c.ptr->c.cmpMula((*(L_DescriptorPtrNodo **)p2)->c.ptr->c);
}

void L_Descriptor3D::crear(const L_Descriptor &d, const L_CoordsCart3D &pos, const L_CoordsCart3D &normal, const L_CoordsCart3D &horz, double distAde, const L_CamaraPinhole &cam)
{
	// Los datos intrinsecos del descriptor x0, y0, sigma0, ang dejan d ser validos
	// Van a ser recalculados cuando se pida el descriptor 3d
	L_CoordsCart3D vert;
	double c, s;
	this->d = d;
	this->d.d3d = NULL;  // Se addFrom al pedir el descriptor
	this->cam = cam;
	this->d.xAde = distAde;
	sigma1 = d.sigma0 * distAde;
	cen.x = pos.x;
	cen.y = pos.y;
	cen.z = pos.z;
	nor.x = normal.x;
	nor.y = normal.y;
	nor.z = normal.z;
	vert = nor.cruz(horz); // Sis de ref local: centro cen, ejes:{normal, horz, vert}

	// El gradiente esta definido en el plano (horz, arr)
	c = cos(d.ang);
	s = sin(d.ang);

	gra.x = horz.x*c + vert.x*s;
	gra.y = horz.y*c + vert.y*s;
	gra.z = horz.z*c + vert.z*s;
}

void L_Descriptor3D::movido_a_fijo(const L_Pose3D_cuat &p, const L_Descriptor3D &orig)
{
	// Recalcula la posicion del descriptor 3d y del descriptor 2d
	*this = orig;  // Copia todo lo necesario
	nor += cen;
	gra += cen;

	cen = p.movido_a_fijo(cen);
	nor = p.movido_a_fijo(nor);
	gra = p.movido_a_fijo(gra);

	nor -= cen;
	gra -= cen;

	d.xAde = cen.x;
	d.ray.tanIzq = cen.y/cen.x;
	d.ray.tanArr = cen.z/cen.x;
	d.d3d = NULL; // Se debe agregar despues
	// Calcular los 4 parametros del descriptor: x0, y0, sigma0, ang
	cam.calcPixel(d.ray, d.x0, d.y0);
	d.sigma0 = sigma1/d.xAde;
	d.ang = atan2(gra.z, -gra.y);
}

void L_Descriptor3D::movido_a_fijo_normalImprecisa(const L_Pose3D_cuat &p, const L_Descriptor3D &orig, double proporcionCorreccion)
{
	// Recalcula la posicion del dscriptor 3d y del descriptor 2d
	L_CoordsCart3D per;
	L_HomogeneousMatrix H;
	L_Pose3D poseRel;
	*this = orig;  // Copia todo lo necesario
	nor += cen;
	gra += cen;

	cen = p.movido_a_fijo(cen);
	nor = p.movido_a_fijo(nor);
	gra = p.movido_a_fijo(gra);

	nor -= cen;
	gra -= cen;
	
	// Girar el sistema de referencia local hacia la camara
	nor = -nor;
	gra = -gra;
	per = nor.cruz(gra);
	H.fija4columnas(nor, gra, per, cen);
	poseRel = H.calcPose3D();
	poseRel.ori.pan *= (1-proporcionCorreccion);
	poseRel.ori.tilt *= (1-proporcionCorreccion);
	H.fijaPose3D(poseRel);
	H.pide4columnas(nor, gra, per, cen);
	nor = -nor;
	gra = -gra;
	
	d.xAde = cen.x;
	d.ray.tanIzq = cen.y/cen.x;
	d.ray.tanArr = cen.z/cen.x;
	d.d3d = NULL; // Se debe agregar despues
	// Calcular los 4 parametros del descriptor: x0, y0, sigma0, ang
	cam.calcPixel(d.ray, d.x0, d.y0);
	d.sigma0 = sigma1/d.xAde;
	d.ang = atan2(gra.z, -gra.y);
}

void L_Descriptor3D::fijo_a_movido(const L_Pose3D_cuat &p, const L_Descriptor3D &orig)
{
	// Recalcula la posicion del descriptor 3d y del descriptor 2d
	*this = orig;  // Copia todo lo necesario
	nor += cen;
	gra += cen;

	cen = p.fijo_a_movido(cen);
	nor = p.fijo_a_movido(nor);
	gra = p.fijo_a_movido(gra);

	nor -= cen;
	gra -= cen;
	// Calcular los 4 parametros del descriptor
	d.ray.tanIzq = cen.y/cen.x;
	d.ray.tanArr = cen.z/cen.x;
	d.xAde = cen.x;
	cam.calcPixel(d.ray, d.x0, d.y0);
	d.sigma0 = sigma1/d.xAde;
	d.ang = atan2(gra.z, -gra.y);
}

void L_Descriptor3DArreglo::generarDescriptorLista(L_DescriptorLista &desL, L_Pose3D_cuat &pCam)
{
	L_Descriptor3D d3d;
	int i;
	for (i=0; i<(int)size(); i++)
	{
		d3d.movido_a_fijo(pCam, (*this)[i]);
		d3d.d.d3d = &operator[](i);
		if (d3d.nor.x > -0.2) // Apunta hacia atras, no agregarlo
			continue;
		desL.push_back(d3d.d);
	}
}

void L_Descriptor3DArreglo::generarDescriptorLista_normalImprecisa(L_DescriptorLista &desL, L_Pose3D_cuat &pCam)
{
	L_Descriptor3D d3d;
	int i;
	for (i=0; i<(int)size(); i++)
	{
		d3d.movido_a_fijo_normalImprecisa(pCam, operator[](i));
		d3d.d.d3d = &operator[](i);
		if (d3d.nor.x > -0.2) // Apunta hacia atras, no agregarlo
			continue;
		desL.push_back(d3d.d);
	}
}

void L_Descriptor3DArreglo::generarDescriptorListaOriginal(L_DescriptorLista &desL)
{
	L_Descriptor3D d3d;
	int i;
	for (i=0; i<(int)size(); i++)
	{
		d3d = (*this)[i];
		d3d.d.d3d = &(*this)[i];
		desL.push_back(d3d.d);
	}
}

void L_Descriptor3DArreglo::genDrawing(L_ShapeArray &lins, L_Pose3D_cuat &pCam)
{
	L_DescriptorLista desL;
	L_DescriptorNodo *desn;
	generarDescriptorLista(desL, pCam);
	int u=0;
	double sigma;
	for (desn=desL.root; desn!=NULL; desn=desn->sig)
	{
		sigma = sqrt(desn->c.sigma0);
		lins.drawArrow((int)(desn->c.x0), (int)(desn->c.y0), sigma, desn->c.ang, (L_uchar)((u*10)%255), (L_uchar)((u*25)%255), (L_uchar)((u*100)%255));
		u++;
	}
}

void L_Descriptor3DArreglo::leerArchivoKb(L_CamaraPinhole &cam, L_Array<L_HomogeneousMatrix> &arrH, std::vector<int> &poses, const char *name, const kbParams &params)
{
	FILE *fp;
	L_Matrix K(3,3);
	L_Pose3D_cuat p3d;
	L_CoordsCart3D ade, izq, arr, tra, cen, nor, horz, posRel;
	L_HomogeneousMatrix Hinv;
	L_Descriptor d;
	int i, j, nt, tmp, numPose;
	bool test = false;

	fp = fopen(name, "r");
	if (fp == NULL)
	{
		printf("No se encuentra el archivo %s\n", name);
		return;
	}
	if (fscanf(fp, "%lf %lf %lf", &K(0,0), &K(0,1), &K(0,2)) != 3)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
	if (fscanf(fp, "%lf %lf %lf", &K(1,0), &K(1,1), &K(1,2)) != 3)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
	if (fscanf(fp, "%lf %lf %lf", &K(2,0), &K(2,1), &K(2,2)) != 3)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}

	cam.fijaMatrizCalibracion(K);

	if (fscanf(fp, "%d", &tmp) != 1)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
	arrH.resize(tmp);

	for (i=0; i<(int)arrH.size(); i++)
	{
		if (fscanf(fp, "%lf %lf %lf %lf", &arrH[i](0,0), &arrH[i](0,1), &arrH[i](0,2), &arrH[i](0,3)) != 4)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		if (fscanf(fp, "%lf %lf %lf %lf", &arrH[i](1,0), &arrH[i](1,1), &arrH[i](1,2), &arrH[i](1,3)) != 4)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		if (fscanf(fp, "%lf %lf %lf %lf", &arrH[i](2,0), &arrH[i](2,1), &arrH[i](2,2), &arrH[i](2,3)) != 4)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		arrH[i].cambiaEjes_transf(arrH[i], params.sis, L_AdeIzqArr);
	}

	if (test)
	{
		K.print("K");
		for (i=0; i<(int)arrH.size(); i++)
			arrH[i].print("H");
	}

	int ene;
	if (fscanf(fp, "%d", &ene) != 1)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
	resize(ene);
	if (fscanf(fp, "%d", &nt) != 1)
		{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
	poses.resize(size());

	d.imorig = NULL;
	d.refGen = NULL;
	d.cam = 0;
	d.nt = nt;
	d.nobj = 0;
	d.d3d = NULL;
	for (i=0; i<(int)size(); i++)
	{
		if (fscanf(fp, "%lf %lf %lf", &cen.x, &cen.y, &cen.z) != 3)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		if (fscanf(fp, "%lf %lf", &d.x0, &d.y0) != 2)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		if (fscanf(fp, "%lf %lf", &d.sigma0, &d.ang) != 2)
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		cen.cambiaEjes(cen, params.sis, L_AdeIzqArr);
		if (params.invertirOri == true)
			d.ang = -d.ang;
		//
		for (j=0; j<nt; j++)
			if (fscanf(fp, "%lf", &d.vector[j]) != 1)
				{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		if (fscanf(fp, "%d", &numPose) != 1) // Se usa esto por mientras
			{printf("L_Descriptor3DArreglo::leerArchivoKb() : formato de archivo incorrecto\n"); fclose(fp); return;}
		poses[i] = numPose;
		// arrH es un sistema de referencia solidario al object cuando el object esta girado
		// sin embargo, ahora el object esta centrado
		Hinv.inversaTrasponiendoDe(arrH[numPose]);
		Hinv.pide4columnas(ade, izq, arr, tra);
		nor = -ade;
		horz = -izq;
		posRel = arrH[numPose].movido_a_fijo(cen);
		operator[](i).crear(d, cen, nor, horz, posRel.x, cam);
	}

	if (test)
	{
		for (i=0; i<(int)size(); i++)
		{
			printf("%g %g %g %ggr -> ", operator[](i).d.x0, operator[](i).d.y0, operator[](i).d.sigma0, operator[](i).d.ang*180/M_PI);
			p3d = arrH[poses[i]].calcPose3D_cuat();
			L_Descriptor3D des;
			des.movido_a_fijo(p3d, (*this)[i]);
			printf("%g %g %g %ggr\n", des.d.x0, des.d.y0, des.d.sigma0, des.d.ang*180/M_PI);
			p3d.imprimeGrados("Pose");
		}
	}

	fclose(fp);
}

void L_Descriptor3DArreglo::calcPose(L_Pose3D_cuat &pose, std::vector<L_Descriptor> &dIm)
{
	L_Matrix uv;
	L_RayoCamara r;
	L_LevenbergMarquardt opt;
	L_Descr3D_datos dat;
	L_Matrix x;

	dat.arreglo = this;
	dat.uv = &uv;

	uv.reallocate((int)dIm.size(), 2);
	for (int i=0; i<uv.li; i++)
	{
		operator[](i).cam.calcRayo(r, dIm[i].x0, dIm[i].y0);
		uv(i,0) = r.tanIzq;
		uv(i,1) = r.tanArr;
	}

	opt.lenVector = 7;
	opt.nIterationsMax = 8;
	x.reallocate(7,1);
	x(0,0) = pose.pos.x;
	x(1,0) = pose.pos.y;
	x(2,0) = pose.pos.z;
	x(3,0) = pose.ori.a;
	x(4,0) = pose.ori.b;
	x(5,0) = pose.ori.c;
	x(6,0) = pose.ori.d;
	opt.minimize_vect((const void *)&dat, x, &errProy);
	pose.pos.x = x(0,0);
	pose.pos.y = x(1,0);
	pose.pos.z = x(2,0);
	pose.ori.a = x(3,0);
	pose.ori.b = x(4,0);
	pose.ori.c = x(5,0);
	pose.ori.d = x(6,0);
}

void L_Descriptor3DArreglo::errProy(const void *data, const L_Matrix &pq, L_Matrix &err)
{
	L_Pose3D_cuat eta;
	L_CoordsCart3D pos;
	L_Descr3D_datos *datos = (L_Descr3D_datos *)data;
	err.reallocate(2*datos->arreglo->size(), 1);
	eta.pos.x = pq(0,0);
	eta.pos.y = pq(1,0);
	eta.pos.z = pq(2,0);
	eta.ori.a = pq(3,0);
	eta.ori.b = pq(4,0);
	eta.ori.c = pq(5,0);
	eta.ori.d = pq(6,0);

	for (int i=0; i<(int)datos->arreglo->size(); i++)
	{
		pos = eta.fijo_a_movido((*datos->arreglo)[i].cen);
		err(2*i+0,0) = datos->uv->operator()(i,0) - pos.y/pos.x;
		err(2*i+1,1) = datos->uv->operator()(i,1) - pos.z/pos.x;
	}
}


void L_Track::createFrom(const L_Descriptor &des, long frame, long idNum, L_uchar R, L_uchar G, L_uchar B, double tiem)
{
	id = idNum;
	d=des;
	cR = R;
	cG = G;
	cB = B;
	x = des.x0;
	y = des.y0;
	t.resize(t.size()+1);
	t[t.size()-1].define(des, frame, (tiem>0) ? tiem : frame);
}

void L_Track::addFrom(const L_Descriptor &des, long frame, double tiem)
{
	d=des;
	t.resize(t.size()+1);
	t.adjustMemDown();
	t[t.size()-1].define(des, frame, (tiem>0) ? tiem : frame);
	x = des.x0;
	y = des.y0;
}

void L_Track::agrega_c2(const L_Descriptor &des, long frame, double tiem)
{
	d=des;
	t.resize(t.size()+1);
	t[t.size()-1].define_c2(des, frame, (tiem>0) ? tiem : frame);
	xc2 = des.x0;
	yc2 = des.y0;
}

void L_Track::addFrom(const L_Descriptor &des, const L_Descriptor &desc2, long frame, double tiem)
{
	d=des;
	t.resize(t.size()+1);
	t[t.size()-1].define(des, frame, (tiem>0) ? tiem : frame);
	t[t.size()-1].define_c2(desc2, frame, (tiem>0) ? tiem : frame);
	x = des.x0;
	y = des.y0;
	xc2 = desc2.x0;
	yc2 = desc2.y0;
}

void L_ConjuntoTracks::procesa1(const L_DescriptorLista &desL, long frame)
{
	L_DescriptorNodo *des;
	int i, iMin;
	double dv, dv1, dv2;
	if (!inicializado()) {printf("L_ConjuntoTracks::procesa() llamado sin init el object\n"); return;}

	// Eliminacion de tracks viejos
	for (i=0; i<(int)arr.size(); i++)
	{
		if (frame - arr[i].t[ arr[i].t.size() - 1 ].frame > frameMax)
		{
			arr[i] = arr[arr.size()-1]; // Se sobrescribe el elemento actual
			arr.resize(arr.size()-1);
			i--; // Para descontar esta iteracion del for
		}
	}

	int arrn=arr.size(); // Si se crean nuevos tracks en esta iteracion no tiene sentido que se le agreguen others descriptores de esta iteracion

	for (des=desL.root; des!=NULL; des=des->sig)
	{
		dv1 = 1e9;
		dv2 = 1e10;
		iMin = -1;
		for (i=0; i<arrn; i++)
		{
			dv = des->c.distEuclidCuadr(arr[i].d);
			if (dv < dv2)
			{
				if (dv < dv1)
				{
					dv2 = dv1;
					dv1 = dv;
					iMin = i;
				}
				else
				{
					dv2 = dv;
				}
			}
		}
		if (dv1 < maxDist && dv1 / dv2 <= maxRatio && fabs(des->c.x0 - arr[iMin].d.x0) < dx &&
			fabs(des->c.y0 - arr[iMin].d.y0) < dy && fabs( L_difang(des->c.ang - arr[iMin].d.ang) ) < dtheta &&
			des->c.sigma0/arr[iMin].d.sigma0 > dsigmaRatio && arr[iMin].d.sigma0/des->c.sigma0 > dsigmaRatio &&
			arr[iMin].t[ arr[iMin].t.size() - 1 ].frame < frame)
		{
			arr[iMin].addFrom(des->c, frame);
		}
		else if ((int)(int)arr.size() < nMax)
		{
			arr.resize(arr.size()+1);
			arr[arr.size()-1].createFrom(des->c, frame, ultId++);
		}
	}
}

void L_ConjuntoTracks::procesa2(const L_DescriptorLista &desL, long frame, L_ImageRGBUchar *im)
{
	L_DescriptorNodo *des;
	int i, iMin;
	double ti, dt, distxy, dv, dv1, dv2, dang, sratio;
	bool bDist, agregado;

	if (!inicializado()) {printf("L_ConjuntoTracks::procesa() llamado sin init el object\n"); return;}

	// Eliminacion de tracks viejos
	ti = L_TIME();
	for (i=0; i<(int)arr.size(); i++)
	{
		dt = ti - arr[i].t[arr[i].t.size()-1].ti;
		//arr[i].x = arr[i].t[arr[i].t.size()-1].x + arr[i].vx*dt;
		//arr[i].y = arr[i].t[arr[i].t.size()-1].y + arr[i].vy*dt;

		if (frame - arr[i].t[ arr[i].t.size() - 1 ].frame > frameMax)
		{
			arr[i] = arr[arr.size()-1]; // Se sobrescribe el elemento actual
			arr.resize(arr.size()-1);
			i--; // Para descontar esta iteracion del for
		}
	}

	int arrn=arr.size(); // Si se crean nuevos tracks en esta iteracion no tiene sentido que se le agreguen others descriptores de esta iteracion

	for (des=desL.root; des!=NULL; des=des->sig)
	{
		dv1 = 1e9;
		dv2 = 1e10;
		iMin = -1;
		agregado = false;
		
		if (arrn > 0)
		{
			// Buscar el mas parecido
			for (i=0; i<arrn; i++)
			{
				distxy = sqrt( (des->c.x0 - arr[i].x)*(des->c.x0 - arr[i].x) + (des->c.y0 - arr[i].y)*(des->c.y0 - arr[i].y) );
				//if ( distxy > 3*dx)
				//	continue;
				if (des->c.vector != NULL)
					dv = des->c.distEuclidCuadr(arr[i].d); // + 0.3*distxy/dx;
				else
					dv = distxy / dx;
				if (dv < dv1)
				{
					dv2 = dv1;
					dv1 = dv;
					iMin = i;
				}
				else if (dv < dv2)
					dv2 = dv;
			}
			if (iMin >= 0)
			{
				bDist = dv1 < maxDist && dv1 / dv2 <= maxRatio;
				dang = fabs( L_difang(des->c.ang - arr[iMin].d.ang) );
				sratio = arr[iMin].d.sigma0 / des->c.sigma0;
				sratio = L_MIN(sratio, 1/sratio);

				if (bDist && fabs(des->c.x0 - arr[iMin].x) < dx &&
					fabs(des->c.y0 - arr[iMin].y) < dy && dang < dtheta &&
					sratio > dsigmaRatio && arr[iMin].t[ arr[iMin].t.size() - 1 ].frame < frame)
				{
					arr[iMin].addFrom(des->c, frame, ti);
					agregado = true;
				}
			}
		}
		if (agregado == false && (int)(int)arr.size() < nMax)
		{
			L_uchar R=255, G=255, B=255;
			arr.resize(arr.size()+1);
			arr[arr.size()-1].t.resize(0);
			if (im!=NULL)
			{
				R = im->pix(int(des->c.x0),int(des->c.y0),0);
				G = im->pix(int(des->c.x0),int(des->c.y0),1);
				B = im->pix(int(des->c.x0),int(des->c.y0),2);
			}
			arr[arr.size()-1].createFrom(des->c, frame, ultId++, R, G, B, ti);
		}
	}
}


void L_ConjuntoTracks::pideDescriptores(std::vector<L_Descriptor> &ret, long frame)
{
	int i, ene;
	// Subir la memoria a ret en una cantidad apropiada
	ene = (int)ret.size();
	ret.resize(ret.size() + arr.size());
	for (i=0; i<(int)arr.size(); i++)
	{
		if (arr[i].t[arr[i].t.size()-1].frame == frame)
			ret[ene+i]=arr[i].d;
	}
}

void L_ConjuntoTracks::procesaFalso(std::vector<L_CoordsCart3D> &p3d, const L_CamaraPinhole &c, L_Pose3D &pose, long frame, int *nVisibles)
{
	int i, j;
	L_Descriptor d;
	double xpix, ypix;
	bool agregado;
	L_CoordsCart3D p2;
	L_RayoCamara r;
	L_HomogeneousMatrix H;
	d.nt = 1;
	H.fijaPose3D_fijo_a_movido(pose); // del sis base al de la camara, transf inversa
	// Eliminacion de tracks viejos
	for (i=0; i<(int)arr.size(); i++)
	{
		if (frame - arr[i].t[ arr[i].t.size() - 1 ].frame > frameMax)
		{
			arr[i] = arr[arr.size()-1]; // Se sobrescribe el elemento actual
			arr.resize(arr.size()-1);
			i--; // Para descontar esta iteracion del for
		}
	}
	for (i=0; i<(int)p3d.size(); i++)
	{
		p2 = H*p3d[i]; // p2 en el sis de la camara
		r.tanIzq = p2.y/p2.x;
		r.tanArr = p2.z/p2.x;
		c.calcPixel(r, xpix, ypix);
		if (p2.x > 0 && xpix > 0 && xpix < c.pideResX() && ypix > 0 && ypix < c.pideResy())
		{
			if (nVisibles != NULL)
				*nVisibles++;
			d.x0 = xpix;
			d.y0 = ypix;
			d.ang = 0;
			d.sigma0 = 1;
			agregado = false;
			for (j=0; j<(int)arr.size(); j++)
			{
				if (fabs(arr[j].d.vector[0] - i) < 0.1)
				{
					arr[j].addFrom(d, frame);
					agregado = true;
					break;
				}
			}
			if (!agregado && (int)arr.size() < nMax)
			{
				arr.resize(arr.size()+1);
				arr[arr.size()-1].createFrom(d, frame, ultId++);
			}

		}
	}
}


void L_ConjuntoTracks::pide_xy_n_frames(int nFrames, int &nTracks, L_Matrix &xy, long frame)
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

void L_ConjuntoTracks::pide_uv_n_frames(int n, int &nTracks, const L_Camara *c, L_Matrix &uv, long frame)
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

void L_FiltroDescriptores::filtrarNodos(L_DescriptorLista &desL)
{
	L_DescriptorNodo *ptr;
	desL.n=0;
	for ( desL.pult=&desL.root ; *desL.pult!=NULL; )
	{
		if ( criterioEliminacion(*desL.pult) )
		{
			ptr=*desL.pult;
			*desL.pult=(*desL.pult)->sig;
			delete ptr;
		}
		else
		{
			desL.pult=&(*desL.pult)->sig;
			desL.n++;
		}
	}
}

void L_FiltroDesPolig::agregaVertice(int x0, int y0)
{
	if (nVertices>=nVerticesMax)
	{
		printf("L_FiltroDesPolig::agregaVertice : mas vertices que la memoria pedida\n");
		return;
	}
	xVer[nVertices]=x0;
	yVer[nVertices]=y0;
	//xCen=(xCen*nVertices+x0)/(nVertices+1);
	//yCen=(yCen*nVertices+y0)/(nVertices+1);
	nVertices++;
}

bool L_FiltroDesPolig::criterioEliminacion(const L_DescriptorNodo *des)
{
	int i1, i2;
	long l;
	for (i1=0; i1<nVertices-1; i1++)
	{
		if (i1<nVertices-2)
			i2=i1+1;
		else
			i2=0;
		l=-(long)L_CrossProduct2D(xVer[i2]-xVer[i1],yVer[i2]-yVer[i1],des->c.x0-xVer[i1],des->c.y0-yVer[i1],double);
		if ( l < 0 )
			return true; // Eliminar
	}
	return false; // No eliminar
}

bool L_FiltroDesPolig::criterioEliminacion(double x, double y)
{
	int i1, i2;
	long l;
	for (i1=0; i1<nVertices-1; i1++)
	{
		if (i1<nVertices-2)
			i2=i1+1;
		else
			i2=0;
		l=-(long)L_CrossProduct2D(xVer[i2]-xVer[i1],yVer[i2]-yVer[i1],x-xVer[i1],y-yVer[i1],double);
		if ( l < 0 )
			return true; // Eliminar
	}
	return false; // No eliminar
}

void L_CalcDistrSigmas::AgregaLista(const L_DescriptorLista &desL)
{
	L_DescriptorNodo *desN;
	for (desN=desL.root; desN!=NULL; desN=desN->sig)
		nPun[desN->c.pir_n+2]++;
}

void L_CalcDistrSigmas::print(bool imprDiv)
{
	int i;
	double sum=0;
	for (i=0; i<50; i++)
		sum+=nPun[i];
	if (imprDiv)
	{
		for (i=0; i<12; i++)
		{
			if (nPun[i]!=0)
				printf("%.1f  ",nPun[i+1]*1.0/nPun[i]);
			else
				printf("### ");
		}
		printf("\n");
	}
	else
	{
		for (i=0; i<12; i++)
			printf("%ld  ",nPun[i]);
		printf("\n");
	}
}

int L_Calce::qcmpRef(const void *calceNodopptr1, const void *calceNodopptr2)
{
	L_CalceNodo *p1 = *(L_CalceNodo **)(calceNodopptr1);
	L_CalceNodo *p2 = *(L_CalceNodo **)(calceNodopptr2);
	return p1->c.dRef->cmpPtr(*p2->c.dRef);
}


int L_Calce::qcmpPosSigmaAng(const void *calceNodopptr1, const void *calceNodopptr2)
{
	L_CalceNodo *p1 = *(L_CalceNodo **)(calceNodopptr1);
	L_CalceNodo *p2 = *(L_CalceNodo **)(calceNodopptr2);
	return p1->c.dRef->cmpPosSigmaAng(*p2->c.dRef);
}
void L_CalceLista::insertaNodoAlInicio_(L_DescriptorLista &desLRef, L_DescriptorLista &desLPru, int xRef, int yRef, int xPru, int yPru, L_ImageGrayDouble *imRef, L_ImageGrayDouble *imPru)
{
	L_Descriptor p1;
	L_Calce cal;
	p1.x0=xRef;
	p1.y0=yRef;
	p1.imorig=imRef;
	p1.imorig_lx=imRef->lx;
	p1.imorig_ly=imRef->ly;
	desLRef.push_front(p1);
	p1.x0=xPru;
	p1.y0=yPru;
	p1.imorig=imPru;
	p1.imorig_lx=imPru->lx;
	p1.imorig_ly=imPru->ly;
	desLPru.push_front(p1);
	cal.dRef = &desLRef.root->c;
	cal.dPru = &desLPru.root->c;
	push_front(cal);
}

void L_CalceLista::genCalces(L_ConjDescrBus &db, L_DescriptorLista &desL, bool restrDeltaAng, double deltaAng, L_List<L_Calce> *cal_pool)
{
	L_DescriptorNodo *ptr;
	std::vector<L_Descriptor *> dest(2);
	std::vector<double> dist(2);
	L_Calce cal;
	bool errTipos = false;

	for (ptr=desL.root; ptr!=NULL; ptr=ptr->sig)
	{
		db.findNearest_anyObject(ptr->c, 2, dest, dist);
		if (dist[1]!=0 && dist[0]<=ptr->c.refGen->calces_maxDist && dist[0]/dist[1]<=ptr->c.refGen->calces_maxRatio)
		{
			if (dest[0]->tipoP != ptr->c.tipoP || dest[0]->tipoD != ptr->c.tipoD)
			{
				errTipos = true;
				continue;
			}
			if (restrDeltaAng)
			{
				double da=ptr->c.ang-dest[0]->ang;
				if (da>M_PI)
					da-=2*M_PI;
				else if (da<-M_PI)
					da+=2*M_PI;
				da=L_ABS(da);
				if (da > deltaAng)
					continue;
			}
			if (ptr->c.intens != dest[0]->intens) // Por si los dos valen zero
			{
				if (ptr->c.signo != dest[0]->signo)
					continue;
				if (ptr->c.intens / dest[0]->intens < ptr->c.refGen->calces_maxDifIntens  ||  dest[0]->intens / ptr->c.intens < ptr->c.refGen->calces_maxDifIntens)
					continue;
			}
			cal.dRef=dest[0];
			cal.dPru = &ptr->c;
			cal.vectDist = dist[0];
			cal.vectDivDist = dist[0]/dist[1];
			push_back(cal, cal_pool);
		}
	}
	if (errTipos)
		printf("L_CalceLista::genCalces() : comparacion de descriptores de distinto tipo\n");
}


#ifdef L_MAX_CALCES_ARR
#error Error: L_MAX_CALCES_ARR ya definido
#endif
#define L_MAX_CALCES_ARR 200
void L_CalceLista::genCalces_nMasCercanos(L_ConjDescrBus &db, L_DescriptorLista &desL, bool restrDeltaAng, double deltaAng, L_List<L_Calce> *cal_pool)
{
	std::vector<L_Descriptor *> dest(L_MAX_CALCES_ARR+1); // Pedir la memoria una sola vez es mas rapido
	std::vector<double> dist(L_MAX_CALCES_ARR+1);
	L_DescriptorNodo *ptr;
	L_Calce cal;
	int i;
	int nMax;
	bool errTipos = false;

	for (ptr=desL.root; ptr!=NULL; ptr=ptr->sig)
	{
		db.findNearest_anyObject(ptr->c, ptr->c.refGen->calces_numCalcesPorDescr+1, dest, dist);
		nMax=L_MIN(L_MAX_CALCES_ARR, ptr->c.refGen->calces_numCalcesPorDescr);
		for (i=0; i<nMax; i++)
		{
			if (dest[i]->tipoP != ptr->c.tipoP || dest[i]->tipoD != ptr->c.tipoD)
			{
				errTipos = true;
				continue;
			}
			if (dist[i+1]!=0 && dist[i]<=ptr->c.refGen->calces_maxDist && dist[i]/dist[nMax]<=ptr->c.refGen->calces_maxRatio)
			{
				if (restrDeltaAng)
				{
					double da=ptr->c.ang-dest[i]->ang;
					if (da>M_PI)
						da-=2*M_PI;
					else if (da<-M_PI)
						da+=2*M_PI;
					da=L_ABS(da);
					if (da > deltaAng)
						continue;
				}
				if (ptr->c.intens != dest[i]->intens) // Por si los dos valen zero
				{
					if (ptr->c.sigma0 != dest[i]->signo)
						continue;
					if (ptr->c.intens / dest[i]->intens < ptr->c.refGen->calces_maxDifIntens  ||  dest[i]->intens / ptr->c.intens < ptr->c.refGen->calces_maxDifIntens)
						continue;
				}
				cal.dRef=dest[i];
				cal.dPru=&ptr->c;
				cal.vectDist = dist[0];
				cal.vectDivDist = dist[0]/dist[1];
				push_back(cal, cal_pool);
			}
			else
				break;
		}
	}
	if (errTipos)
		printf("L_CalceLista::genCalces() : comparacion de descriptores de distinto tipo\n");
}
#undef L_MAX_CALCES_ARR

double L_CalceLista::genCalcesVerif(L_BusqSecDescr &sec, L_ConjDescrBus &db, L_DescriptorLista &desL, double distmax, double ratio)
{
	int ntot=0, ncorr=0;
	L_DescriptorNodo *ptr;
	std::vector<L_Descriptor *> dest(2);
	std::vector<double> dist(2);
	L_Calce cal;

	for (ptr=desL.root; ptr!=NULL; ptr=ptr->sig)
	{
		db.findNearest_anyObject(ptr->c, 2, dest, dist);
		if (dist[1]!=0 && dist[0]<=distmax && dist[0]/dist[1]<=ratio)
		{
			ntot++;
			if ( sec.verifBusqCorrecta(ptr->c, *dest[0], *dest[1]) == L_BusqSecDescr::sin_error)
				ncorr++;
			cal.dRef=dest[0];
			cal.dPru=&ptr->c;
			cal.vectDist = dist[0];
			cal.vectDivDist = dist[0]/dist[1];
			push_back(cal);
		}
	}
	return ncorr/(double)ntot;
}

void L_CalceLista::parear(L_DescriptorLista &desL1, L_DescriptorLista &desL2)
{
	L_DescriptorNodo *d1, *d2;
	L_Calce c;

	if (desL1.size() != desL2.size())
	{
		printf("L_CalceLista::parear() : error en length de argumentos\n");
		return;
	}

	for (d1=desL1.root,d2=desL2.root; d1!=NULL&&d2!=NULL; d1=d1->sig,d2=d2->sig)
	{
		c.dRef=&d1->c;
		c.dPru=&d2->c;
		push_back(c);
	}
}

bool L_CalceLista::dibuja(const char *nomarch) //! Solo puede dibujar si hay calces
{
	L_ImageRGBUchar im;
	L_ShapeArray lins;
	L_Array<const L_ImageGrayDouble *> arr(2);
	if (root==NULL)
		return false;
	arr[0] = root->c.dPru->imorig;
	arr[1] = root->c.dRef->imorig;
	genLineas(lins);
	im.genDrawingMatches(arr, lins);
#if __COMPAT_ATLIMAGE__
	return im.writeImageATL(nomarch)!=0;
#else
	return im.writeBMP(nomarch)!=0;
#endif
}

bool L_CalceLista::dibuja(L_ImageRGBUchar &im) //! Solo puede dibujar si hay calces
{
	if (root==NULL)
		return false;
	L_ShapeArray lins;
	L_Array<const L_ImageGrayDouble *> arr(2);
	arr[0] = root->c.dPru->imorig;
	arr[1] = root->c.dRef->imorig;
	genLineas(lins);
	im.genDrawingMatches(arr, lins);
	return true;
}

void L_CalceLista::dibuja(L_ImageRGBUchar &im, L_ImageGrayDouble &imRef, L_ImageGrayDouble &imPru) //! Solo puede dibujar si hay calces
{
	L_Array<const L_ImageGrayDouble *> arr(2);
	L_ShapeArray lins;
	if (root==NULL)
	{
		arr[0] = &imPru;
		arr[1] = &imRef;
	}
	else
	{
		arr[0] = root->c.dPru->imorig;
		arr[1] = root->c.dRef->imorig;
	}
	genLineas(lins);
	im.genDrawingMatches(arr, lins);
}

void L_SIFTdireccLista::calcDirecciones(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin)
{
	if (orig)
		L_SIFTdireccLista::calcDirecciones_Lowe(xc, yc, grad1, grad2, gradPolares, s, pin);
	else
		L_SIFTdireccLista::calcDirecciones_PLZ(xc, yc, grad1, grad2, gradPolares, s, pin);
}

void L_SIFTdireccLista::calcDirecciones_PLZ(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin)
{
#ifdef L_BASIC_DEBUG
	if (nd<=0)
		L_hard_shutdown("L_SIFTdireccLista::calcDirecciones_PLZ : non-positive size allocation");
#endif
	std::vector<double> bins(nd), bins2(nd);
	int i, j;
	int ind;
	int s3;
	double grmod, grang;
	L_DireccNodo dnodo;
	double x0, x1, x2, y0, y1, y2, ang;
	double Fac = 1/(2*M_PI)*nd;

	for (i=0; i<nd; i++)
		bins[i]=0;

	s3=L_GaussianHalfWidth(s*facS);

	if (gaus.size() == 0)
		L_calcGausHNorm(gaus,s*facS,2*s3+1);
	if (xc-s3<0 || xc+s3>=grad1.lx || yc-s3<0 || yc+s3>=grad1.ly) // punto muy cerca del borde
		return;

	if (gradPolares) // Si el gradiente ya esta descompuesto en coordenadas polares, aprovechar eso
	{
		for (j=yc-s3; j<=yc+s3; j++)
		{
			for (i=xc-s3; i<=xc+s3; i++)
			{
				grmod=grad1.pix(i,j);
				grang=grad2.pix(i,j);
				ind=L_FLOOR1000(grang*Fac) % nd;
				if (ind<0) ind+=nd;
				if (ind>=nd) ind-=nd;
				bins[ind]+=grmod*gaus[i-xc+s3]*gaus[j-yc+s3];
			}
		}
	}
	else
	{
		for (j=yc-s3; j<=yc+s3; j++)
		{
			for (i=xc-s3; i<=xc+s3; i++)
			{
				grmod=sqrt(grad1.pix(i,j)*grad1.pix(i,j)+grad2.pix(i,j)*grad2.pix(i,j));
				grang=atan2(grad2.pix(i,j),grad1.pix(i,j));
				ind=L_FLOOR1000(grang*Fac) % nd;
				if (ind<0) ind+=nd;
				if (ind>=nd) ind-=nd;
				bins[ind]+=grmod*gaus[i-xc+s3]*gaus[j-yc+s3];
			}
		}
	}

	for (i=0; i<nd; i++)
		bins2[i] = (0.25*bins[(i-1+nd)%nd] + 0.5*bins[i] + 0.25*bins[(i+1)%nd]);
	for (i=0; i<nd; i++)
		bins[i] = (0.25*bins2[(i-1+nd)%nd] + 0.5*bins2[i] + 0.25*bins2[(i+1)%nd]);

	// Busqueda de la direccion principal
	ind=0;
	for (i=1; i<nd; i++)
	{
		if (bins[i]>=bins[ind])
			ind=i;
	}
	// Busqueda de direcciones dominantes, interpolacion del maxElement
	for (i=0; i<nd; i++)
	{
		y0 = bins[(i-1+nd)%nd];
		y1 = bins[i];
		y2 = bins[(i+1)%nd];

		if (y1<=y0 || y1<=y2) // no es máximo
			continue;
		if (y1<porc*bins[ind])
			continue;

		// calcular dirección principal del máximo mediante interpolación parabólica
		x0=(0.5+i-1)*2*M_PI/nd;  // posición central del sector angular correspondiente a k-1
		x1=(0.5+i)*2*M_PI/nd;  // posición central del sector angular correspondiente a k
		x2=(0.5+i+1)*2*M_PI/nd;  // posición central del sector angular correspondiente a k+1

		ang=-x0*y2+x0*y1-y1*x2+x1*y2-y0*x1+y0*x2; // denominador

		if (ang==0) // DEBUG : no debería ocurrir
			continue;

		ang=0.5*(y2*x1*x1-x1*x1*y0+y0*x2*x2+y1*x0*x0-x0*x0*y2-y1*x2*x2)/ang;

		if (ang<x0 || ang>x2) // DEBUG : no debería ocurrir
			continue;

		if (ang>=2*M_PI) ang-=2*M_PI;
		if (ang<0) ang+=2*M_PI;

		dnodo = *pin;
		dnodo.c.ang=ang;
		dnodo.sig=NULL;
		lis.push_back(dnodo.c);
	}
	return;
}

void L_SIFTdireccLista::calcDirecciones_Lowe(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin)
{
//Float FindOri(Image grad, Image ori, int row, int col)
	int i, r, c, rows, cols, radius, bin;
	double distsq, gval, weight, angle;
	L_Direcc dnodo;
	double OriSigma = facS * s; // valía 3

	std::vector<double> hist(nd);

	rows = grad1.ly;
	cols = grad1.lx;

	for (i = 0; i < nd; i++)
	hist[i] = 0.0;

	// Look at pixels within 3 sigma around the point and put their
	// Gaussian weighted values in the histogram.
	radius = (int) (OriSigma * 3.0);
	if (xc <= radius || xc >= grad1.lx-radius || yc <= radius || yc >=grad2.ly-radius)
		return;
	if (gradPolares)
	{
		for (r = yc - radius; r <= yc + radius; r++)
			for (c = xc - radius; c <= xc + radius; c++)
				// Do not use last yc or column, which are not valid.
				if (r >= 0 && c >= 0 && r < rows - 2 && c < cols - 2) {
					gval = grad1.pix(c,r);
					distsq = (r - yc) * (r - yc) + (c - xc) * (c - xc);
					if (gval > 0.0  &&  distsq < radius * radius + 0.5) {
						weight = exp(- distsq / (2.0 * OriSigma * OriSigma));
						// Ori is in range of -PI to PI.
						angle = grad2.pix(c,r);
						bin = (int) (nd * (angle + M_PI + 0.001) / (2.0 * M_PI));
						if ( ! (bin >= 0 && bin <= nd)  )
						{
							printf("Error\n");
							return;
						}
						bin = L_MIN(bin, nd - 1);
						hist[bin] += weight * gval;
					}
				}
	}
	else
	{
		for (r = yc - radius; r <= yc + radius; r++)
			for (c = xc - radius; c <= xc + radius; c++)
				// Do not use last yc or column, which are not valid.
				if (r >= 0 && c >= 0 && r < rows - 2 && c < cols - 2) {
					gval = sqrt(grad1.pix(c,r)*grad1.pix(c,r) + grad2.pix(c,r)*grad2.pix(c,r));
					distsq = (r - yc) * (r - yc) + (c - xc) * (c - xc);
					if (gval > 0.0  &&  distsq < radius * radius + 0.5) {
						weight = exp(- distsq / (2.0 * OriSigma * OriSigma));
						// Ori is in range of -PI to PI.
						if (grad2.pix(c,r)==0 && grad1.pix(c,r)==0)
							angle = 0; // Correccion PLZ
						else
							angle = atan2(grad2.pix(c,r), grad1.pix(c,r));
						bin = (int) (nd * (angle + M_PI + 0.001) / (2.0 * M_PI));
						if ( ! (bin >= 0 && bin <= nd)  )
						{
							printf("Error\n");
							return;
						}
						bin = L_MIN(bin, nd - 1);
						hist[bin] += weight * gval;
					}
				}
	}
	// Apply smoothing twice.
	SmoothHistogram(hist, nd);
	SmoothHistogram(hist, nd);

	dnodo=pin->c;
	dnodo.ang=FindOriPeaks(hist, nd); // Se multiplica por -1 por sistema de referencia usado
	lis.push_back(dnodo);
	return;
}

void L_SIFTdireccLista::SmoothHistogram(std::vector<double> &hist, int bins)
{
   int i;
   double prev, temp;

   prev = hist[bins - 1];
   for (i = 0; i < bins; i++) {
      temp = hist[i];
      hist[i] = 0.25 * prev + 0.5 * hist[i] + 0.25 *
         hist[(i + 1 == bins) ? 0 : i + 1];
      prev = temp;
   }
}


double L_SIFTdireccLista::FindOriPeaks(std::vector<double> &hist, int bins)
{
   int i, maxloc = 0;
   double maxval = 0.0;

   /* Find peak in histogram. */
   for (i = 0; i < bins; i++)
      if (hist[i] > maxval) {
         maxval = hist[i];
         maxloc = i;
      }
      /* Set angle in range -PI to PI. */
   return (2.0 * M_PI * (maxloc + 0.5) / bins - M_PI);
}

bool L_SIFTCalculador::calcDescriptor(L_ImageGrayDouble &grad1, L_ImageGrayDouble &grad2, bool gradPolares) // Usa ->x e ->y del descriptor
{
	nt = nx*ny*nd;
	if (orig)
		return calcDescriptorGenerico_Lowe(grad1, grad2, gradPolares);
	else
	{
		const L_SIFTGenFnDescriptorNodo *fnPtr;
		fnPtr=L_SIFTGenFnArr.buscaFn(lx, ly, nd, nx, ny);
		if (fnPtr==NULL)
		{
			return calcDescriptorGenerico_norot(grad1, grad2, gradPolares);
		}
		else
		{
			return (this->*(fnPtr->f))(grad1, grad2, gradPolares);
		}
	}
}

bool L_SIFTCalculador::calcDescriptorGenerico_Lowe(L_ImageGrayDouble &grad1, L_ImageGrayDouble &grad2, bool gradPolares)
{
//void MakeKeypointSample(Keypoint key, Image grad, Image ori, int pool)
	int i;
	/* Produce sample vector. */
	KeySample(vector, grad1, grad2, (int)pir_y, (int)pir_x, gradPolares);

	/* Normalize vector.  This should provide illumination invariance for
	planar lambertian surfaces (except for saturation effects).   */
	NormalizeVec(vector, nt);

	/* Now that illumination normalization has been done, threshold elements
	of index vector to decrease emphasis on large gradient magnitudes. */
	if (nd != 2) {
		for (i = 0; i < nt; i++)
			if (vector[i] > 0.2)//MaxIndexVal)
				vector[i] = 0.2;//MaxIndexVal;
		/* Again normalize vector.  This slightly improves nearest-neighbor
		by increasing relative distance for vectors with few features.  This
		is also useful to implement a distance threshold and to allow
		conversion to integer format. */
		NormalizeVec(vector, nt);
	}

	/* Copy double vector to key vec. */
	return true;
}
void L_SIFTCalculador::NormalizeVec(double *vec, int len)
{
	int i;
	double val, fac, sqlen = 0.0;

	for (i = 0; i < len; i++) {
		val = vec[i];
		sqlen += val * val;
	}
	fac = 1.0 / sqrt(sqlen);
	for (i = 0; i < len; i++)
		vec[i] *= fac;
}
void L_SIFTCalculador::KeySample(double *fvec, L_ImageGrayDouble & grad, L_ImageGrayDouble & ori, int row, int col, bool gradPolares)
{
	int i, j, k, radius, v;
	double sine, cosine, rrot, crot, rloc, cloc;
	double ***index=L_new3d<double>(nx,ny,nd);//[IndexSize][IndexSize][ang];

	/* Initialize index array. */
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			for (k = 0; k < nd; k++)
				index[i][j][k] = 0.0;

	sine = sin(ang); // ang está medido en contra del eje y
	cosine = cos(ang); // ang está medido en contra del eje y

	/* Radius of index sample region in pixels = sqrt(2)/2 * square width. */
	radius = (int) (L_MAX(nx,ny) * (double) L_MAX(lx,ly) / 1.414 + 0.5);

	/* Examine all points from the gradient image that could lie within the
	  index square. */
	for (i = -radius; i <= radius; i++)
		for (j = -radius; j <= radius; j++) {

		/* Rotate sample offset to make it relative to key orientation.
		Uses (row,col) instead of (x,y) coords. */
		//rrot = cosine * i + sine * j;
		//crot = - sine * i + cosine * j;

		crot = cosine * j - sine * i;
		rrot = cosine * i + sine * j;

		/* Compute location of sample in terms of real-valued index
		coordinates.  Subtract 0.5 so that rloc of 1.0 means to put full
		weight on index[1] (eg, when rrot is 0 and IndexSize is 3). */
		rloc = rrot / ly + ny / 2.0 - 0.5;
		cloc = crot / lx + nx / 2.0 - 0.5;

		/* Test whether this sample falls within boundary of the index. */
		if (rloc > -1.0 && rloc < (double) ny  &&
		cloc > -1.0 && cloc < (double) nx)
			AddSample(index, grad, ori, row + i, col + j, rloc, cloc, gradPolares);
	}

	/* Copy the index into a vector attached to this keypoint. */
	v = 0;
	for (i = 0; i < nx; i++)
	 for (j = 0; j < ny; j++)
		for (k = 0; k < nd; k++)
			fvec[v++] = index[i][j][k];
	L_delete3d<double>(index);
}


/* Given a sample from the image gradient, place it in the index array.
*/
void L_SIFTCalculador::AddSample(double ***index, L_ImageGrayDouble & grad, L_ImageGrayDouble & orim, int r, int c, double rloc, double cloc, bool gradPolares)
{
	double mag, ori;
	int IgnoreGradSign=false;

	/* Clip at image boundaries (note last row and col of grad is empty). */
	if (r < 0  ||  r >= (grad.ly - 1) || c < 0  ||  c >= (grad.lx - 1))
		return;

	/* If gradient magnitude is below threshold, then no need to index. */
	if (gradPolares)
		mag = grad.pix(c,r);
	else
		mag = sqrt(grad.pix(c,r)*grad.pix(c,r) + orim.pix(c,r)*orim.pix(c,r));

	/* Subtract keypoint orientation to give ori relative to keypoint. */
	if (gradPolares)
		ori = orim.pix(c,r) - ang;
	else
	{
		if (orim.pix(c,r)==0 && grad.pix(c,r)==0)
			ori=-ang;
		else
			ori = atan2(orim.pix(c,r),grad.pix(c,r)) - ang;
	}

	/* Put orientation in range [0, 2*PI].  If sign of gradient is to
	be ignored, then put in range [0, PI]. */
	if (IgnoreGradSign) {
		while (ori > M_PI)
			ori -= M_PI;
		while (ori < 0.0)
			ori += M_PI;
	} else {
	   while (ori > 2*M_PI)
			ori -= 2*M_PI;
	   while (ori < 0.0)
			ori += 2*M_PI;
	}
	if (nd == 2)
		PlaceInIndexVec(index, mag, ori, rloc, cloc);
	else
		PlaceInIndex(index, mag, ori, rloc, cloc);
}


/* Increment the appropriate locations in the index to incorporate
   this image sample.  The location of the sample in the index is (rloc,cloc).
*/
void L_SIFTCalculador::PlaceInIndex(double ***index, double mag, double ori, double rloc, double cloc)
{
	int r, c, o_r, ri, ci, oi, rindex, cindex, oindex;
	double oval, rfrac, cfrac, ofrac, rweight, cweight, oweight;
	double *ivec;
	bool IgnoreGradSign=false;

	oval = nd * ori / (IgnoreGradSign ? M_PI : 2*M_PI);

	ri = (int)(  (rloc >= 0.0) ? rloc : rloc - 1.0  );  /* Round down to next integer. */
	ci = (int)(  (cloc >= 0.0) ? cloc : cloc - 1.0  );
	oi = (int)(  (oval >= 0.0) ? oval : oval - 1.0  );
	rfrac = rloc - ri;         /* Fractional part of location. */
	cfrac = cloc - ci;
	ofrac = oval - oi;
	if ( ! (ri >= -1  &&  ri < ny  &&  oi >= 0  &&  oi <= nd  && rfrac >= 0.0  &&  rfrac <= 1.0)  )
	{
		printf("Error\n");
		return;
	};

	/* Put appropriate fraction in each of 8 buckets around this point
	  in the (row,col,ori) dimensions.  This loop is written for
	  efficiency, as it is the inner loop of key sampling. */
	for (r = 0; r < 2; r++) {
		rindex = ri + r;
		if (rindex >=0 && rindex < ny) {
			rweight = mag * ((r == 0) ? 1.0 - rfrac : rfrac);

			for (c = 0; c < 2; c++) {
				cindex = ci + c;
				if (cindex >=0 && cindex < nx) {
					cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
					ivec = index[rindex][cindex];

					for (o_r = 0; o_r < 2; o_r++) {
						oindex = oi + o_r;
						if (oindex >= nd)  /* Orientation wraps around at PI. */
							oindex = 0;
						oweight = cweight * ((o_r == 0) ? 1.0 - ofrac : ofrac);
						ivec[oindex] += oweight;
					}
				}
			}
		}
	}
}


/* Like PlaceInIndex, but this sums the (u,v) coordinates of each gradient
   vector rather than using orientation planes.
*/
void L_SIFTCalculador::PlaceInIndexVec(double ***index, double mag, double ori, double rloc, double cloc)
{
	int r, c, ri, ci, rindex, cindex;
	double xgrad, ygrad, rfrac, cfrac, rweight, cweight;
	double *ivec;

	/* Create x,y components of the gradient vector.  In the future, this
	  can be made more efficient by using original gradient values. */
	xgrad = mag * cos(ori);
	ygrad = mag * sin(ori);

	ri = (int)(  (rloc >= 0.0) ? rloc : rloc - 1.0  );  /* Round down to next integer. */
	ci = (int)(  (cloc >= 0.0) ? cloc : cloc - 1.0  );
	rfrac = rloc - ri;         /* Fractional part of location. */
	cfrac = cloc - ci;
	if ( ! (ri >= -1  &&  ri < ny  &&  rfrac >= 0.0  &&  rfrac <= 1.0)  )
	{
		printf("Error\n");
		return;
	};

	/* Put appropriate fraction in each of 4 buckets around this point
	  in the (row,col) dimensions.  This loop is written for
	  efficiency, as it is the inner loop of key sampling. */
	for (r = 0; r < 2; r++) {
		rindex = ri + r;
		if (rindex >=0 && rindex < ny) {
			rweight = ((r == 0) ? 1.0 - rfrac : rfrac);

			for (c = 0; c < 2; c++) {
				cindex = ci + c;
				if (cindex >=0 && cindex < nx) {
					cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
					ivec = index[rindex][cindex];
					ivec[0] += cweight * xgrad;
					ivec[1] += cweight * ygrad;
				}
			}
		}
	}
}







bool L_SIFTCalculador::calcDescriptorGenerico_norot(L_ImageGrayDouble &grad1, L_ImageGrayDouble &grad2, bool gradPolares)
{
	double cosTheta, senTheta;
	double u, v; // posición relativa a la ventana de analisis
	int x, y; // posición absoluta en la imagen
	int nxt, nyt;
	double resang;
	L_LinearQuantizer c_i, c_j, c_k;

	long i[2], j[2], k[2]; // indices
	double du[2], dv[2], dang[2]; // distancias del punto (u,v,ang) actual a los puntos centrales de los bins
	double ucen[2], vcen[2], angcen[2]; // puntos centrales de los bins
	int cont;

	double factorUV;
	double grang;
	double grmod;
	double angLocal;
	double xc, yc;
	int xIni, xFin, yIni, yFin;
	int xrnd, yrnd;

	//matriz de rotacion para "eje v" invertido
	cosTheta=cos(this->ang);
	senTheta=sin(this->ang);
	xc=this->pir_x-0.5; // Para pasar de this->x en la imagen normal a xc en el gradiente hacia adelante
	yc=this->pir_y-0.5; // Para pasar de this->x en la imagen normal a xc en el gradiente hacia adelante

	//definiciones de uso interno: para rotacion
	#if defined(calcU) || defined(calcV) || defined(u_) || defined(v_)
		#error #define conflicts
	#endif
	#define calcU(X,Y) L_ROUND( ((X)-xc)*cosTheta - ((Y)-yc)*senTheta ) // +0.5 para que quede redondeado
	#define calcV(X,Y) L_ROUND( ((X)-xc)*senTheta + ((Y)-yc)*cosTheta )
	#define u_ (L_FLOOR(u))
	#define v_ (L_FLOOR(v))

	nxt=lx*nx/2; // dg_nx DEBE ser par.
	nyt=ly*ny/2; // dg_ny DEBE ser par.
	resang=2*M_PI/nd;  // resolucion angular de los bines del descriptor SIFT

	xIni=L_FLOOR(xc-1.4142*nxt + 0.5);
	if (xIni<0)
		return false;
	xFin=L_FLOOR(xc+1.4142*nxt + 0.5);
	if (xFin>=grad1.lx)
		return false;
	yIni=L_FLOOR(yc-1.4142*nyt + 0.5);
	if (yIni<0)
		return false;
	yFin=L_FLOOR(yc+1.4142*nyt + 0.5);
	if (yFin>=grad1.ly)
		return false;

	for (cont=0; cont<nt; cont++)
		vector[cont]=0;

	if (gdx_h.size() == 0)
		L_calcGausHNorm(gdx_h,nxt,2*nxt); //gdx_h[i]: i debe ir entre (-nxt)+nxt y (nxt-1)+nxt
	if (gdy_h.size() == 0)
		L_calcGausHNorm(gdy_h,nyt,2*nyt); //gdx_h[i]: i debe ir entre (-nyt)+nyt y (nyt-1)+nyt

	c_i.setFromCellRange(-nxt,-nxt+lx,0);
	c_j.setFromCellRange(-nyt,-nyt+ly,0);
	c_k.setFromCellRange(-resang/2,resang/2,0);

	for (y=yIni; y<=yFin; y++)
	{
		for (x=xIni; x<=xFin; x++)
		{
			u=calcU(x,y);
			v=calcV(x,y);

			if (u<-nxt || u>=nxt || v<-nyt || v>=nyt)
				continue;

			i[0]=L_LinearQuantizer_nClusterLeft(c_i,u);
			i[1]=i[0]+1;

			ucen[0]=L_LinearQuantizer_rCenCluster(c_i,i[0]);
			ucen[1]=ucen[0]+lx;

			du[0]=1-L_ABS(u-ucen[0])/lx; //ERROR GRAVE: no tenia denominador
			du[1]=1-L_ABS(u-ucen[1])/lx;

			j[0]=L_LinearQuantizer_nClusterLeft(c_j,v);
			j[1]=j[0]+1;

			vcen[0]=L_LinearQuantizer_rCenCluster(c_j,j[0]);
			vcen[1]=vcen[0]+ly;

			dv[0]=1-L_ABS(v-vcen[0])/ly;
			dv[1]=1-L_ABS(v-vcen[1])/ly;

			xrnd = L_ROUND(x);
			yrnd = L_ROUND(y);

			if (gradPolares)
			{
				grmod=grad1.pix(xrnd,yrnd);
				grang=grad2.pix(xrnd,yrnd);
			}
			else
			{
				grmod=sqrt(grad1.pix(xrnd,yrnd)*grad1.pix(xrnd,yrnd)+grad2.pix(xrnd,yrnd)*grad2.pix(xrnd,yrnd));
				if (grad2.pix(xrnd,yrnd)==0 && grad1.pix(xrnd,yrnd)==0)
					grang=0;
				else
					grang=atan2(-grad2.pix(xrnd,yrnd),grad1.pix(xrnd,yrnd));
			}
			angLocal=grang-this->ang;

			k[0]=L_LinearQuantizer_nClusterLeft(c_k,angLocal);
			k[1]=k[0]+1;  // k[1]*nx*ny

			angcen[0]=L_LinearQuantizer_rCenCluster(c_k,k[0]);
			angcen[1]=angcen[0]+resang;

			L_SETABS(dang[0], (angLocal-angcen[0])/resang);
			L_SETABS(dang[1], (angLocal-angcen[1])/resang);

			/*
			casos:
				para i1, i2:
				* i1<0 : todo a i2
				* (i1, i2) valido : se reparten entre i1 e i2
				* i2>=dg_nx : todo a i1
			Pasa lo mismo con (j1,j2).

			Con (k1,k2) el k es "periodico"
			*/

			if (i[0]==-1) i[0]++;
			if (i[1]>=nx) i[1]--;
			if (j[0]==-1) j[0]++;
			if (j[1]>=ny) j[1]--;
			while (k[0]<0) k[0]+=nd;
			while (k[0]>=nd) k[0]-=nd;
			while (k[1]<0) k[1]+=nd;
			while (k[1]>=nd) k[1]-=nd;
			factorUV=gdx_h[u_+nxt]*gdy_h[v_+nyt]*grmod;

			vector[i[0]+j[0]*nx+k[0]*nx*ny ]+=(double)(
				du[0]*dv[0]*dang[0]*factorUV );

			vector[i[1]+j[0]*nx+k[0]*nx*ny ]+=(double)(
				du[1]*dv[0]*dang[0]*factorUV );

			vector[i[0]+j[1]*nx+k[0]*nx*ny ]+=(double)(
				du[0]*dv[1]*dang[0]*factorUV );

			vector[i[1]+j[1]*nx+k[0]*nx*ny ]+=(double)(
				du[1]*dv[1]*dang[0]*factorUV );

			vector[i[0]+j[0]*nx+k[1]*nx*ny ]+=(double)(
				du[0]*dv[0]*dang[1]*factorUV );

			vector[i[1]+j[0]*nx+k[1]*nx*ny ]+=(double)(
				du[1]*dv[0]*dang[1]*factorUV );

			vector[i[0]+j[1]*nx+k[1]*nx*ny ]+=(double)(
				du[0]*dv[1]*dang[1]*factorUV );

			vector[i[1]+j[1]*nx+k[1]*nx*ny ]+=(double)(
				du[1]*dv[1]*dang[1]*factorUV );

		}
	}

	#undef calcU
	#undef calcV
	#undef u_
	#undef v_

	normalizaHist();
	for (cont=0; cont<nt; cont++)
	{
		if (vector[cont]>0.2)
			vector[cont]=0.2;
		if (vector[cont]<-0.2)
			vector[cont]=-0.2;
	}
	normalizaHist();
	return true;
}

//////////////////////////
//// kd-tree recursivo

// Clases derivadas para kd-tree con busqueda BBF
class L_DescriptorRefOrd
{
public:
	L_DescriptorNodo *nodo; // Una arreglo L_DescriptorRefOrd[] se sort con qsort segun num, no contiene la memoria del object
	double num;
	bool asign; // Indica si este nodo ya ha sido absorbido por el kd-tree
	L_DescriptorRefOrd() {nodo=NULL; num=0; asign=false;}
	static bool sort(L_DescriptorRefOrd *arr, int n, int nc);
	bool distintoA(L_DescriptorRefOrd &other) {return nodo->c.distEuclidCuadr(other.nodo->c) > 0;}
	static int cmp(const void *nodo1, const void *nodo2);
	static int cmpVect(const void *nodo1, const void *nodo2);
};

bool L_DescriptorRefOrd::sort(L_DescriptorRefOrd *arr, int n, int nc)
{
	int i;
	// Definir dimension a comparar en los nodos
	for (i=0; i<n; i++)
		arr[i].num=arr[i].nodo->c.vector[nc];
	// Ordenar nodos
	qsort(arr, n, sizeof(L_DescriptorRefOrd), L_DescriptorRefOrd::cmp);
	return true;
}

#if defined L_DesRefOrdNnum
	#error #define conflicts
#endif
#define L_DesRefOrdNnum(x) ( ((L_DescriptorRefOrd *)(x))->num )
int L_DescriptorRefOrd::cmp(const void *nodo1, const void *nodo2)
{
	return (L_DesRefOrdNnum(nodo1)>L_DesRefOrdNnum(nodo2)) - (L_DesRefOrdNnum(nodo1)<L_DesRefOrdNnum(nodo2));
}
#undef L_DesRefOrdNnum

int L_DescriptorRefOrd::cmpVect(const void *nodo1, const void *nodo2)
{
	return ((L_DescriptorRefOrd *)nodo1)->nodo->c.cmpMula(((L_DescriptorRefOrd *)nodo2)->nodo->c);
	//return (*(L_DescriptorNodo **)nodo1)->cmpMula(**(L_DescriptorNodo **)nodo2);
}

bool L_KdNodo::creaRec(void *LDRF_arr, int n, L_KdTreeBBF *kdt)
{
	L_DescriptorRefOrd *arr = (L_DescriptorRefOrd *)LDRF_arr;
	double maxvar=-1;
	double sum, suma2, var, d, dmin;
	int imax=0, imed;
	int i, j;

	// Un L_KdNodo tiene 3 punteros: izq, der y descr
	// Se puede llenar el par (izq,der) o bien (descr) pero no todos a la vez

	// Voy a ser hoja -> lleno descr
	if (n==1)
	{
		arr[0].nodo->sig = NULL; // Ya no esta en la lista, esta en el tree
		descr = arr[0].nodo;
		arr[0].asign = true;
		return true;
	}

	// No es hoja -> Es nodo de division -> lleno izq y der

	// buscar i de maxima variance => imax
	for (i=0; i<arr[0].nodo->c.nt; i++)
	{
		sum=0;
		suma2=0;
		for(j=0; j<n; j++)
		{
			sum+=arr[j].nodo->c.vector[i];
			suma2+=arr[j].nodo->c.vector[i]*arr[j].nodo->c.vector[i];
		}
		var=suma2-sum*sum/n;
		if (var>maxvar)
		{
			maxvar=var;
			imax=i;
		}
	}
	// imax es la componente de mayor variance
	// ordenar
	L_DescriptorRefOrd::sort(arr, n, imax); // Ordenar arr según la componente imax del descriptor

	// buscar la mediana (en este caso la transicion de valores mas cercana a la mitad del arreglo)
	// -> arr se divide en [0,...,imed-1] y [imed,...,n-1]
	dmin = 1.0e10;
	imed = -1;
	for (i=1; i<n; i++)
	{
		d = n/2 - i;
		if (d<0)
			d = -d;
		if (d < dmin && arr[i-1].num != arr[i].num)
		{
			dmin = d;
			imed = i;
		}
	}
	if (imed == -1) // todos los elementos son iguales y no hay mediana
		return false;

	// Crear el punto de division
	dim=imax;
	umbral=(arr[imed-1].nodo->c.vector[dim] + arr[imed].nodo->c.vector[dim]) / 2.0; // divide [0,...,imed-1] y [imed,...,n-1]
	// Analizar el lado izquierdo [0,...,imed-1]
	izq = new L_KdNodo();
	if (izq->creaRec(&arr[0], imed, kdt) == false)
		return false;
	
	// Analizar el lado derecho [imed,...,n-1]
	der = new L_KdNodo();
	if (der->creaRec(&arr[imed], n-imed, kdt) == false)
		return false;

	return true;
}

//////////////////////////
//// kd-tree no recursivo



int L_KdNodo_cmpRepetidos(const void *kdNodo_1, const void *kdNodo_2)
{
	L_KdNodo *p1 = (L_KdNodo *)kdNodo_1;
	L_KdNodo *p2 = (L_KdNodo *)kdNodo_2;
	int i;
	for (i=0; i<p1->descr->c.nt; i++)
	{
		if (p1->descr->c.vector[i] > p2->descr->c.vector[i])
			return 1;
		else if (p1->descr->c.vector[i] < p2->descr->c.vector[i])
			return -1;
	}
	return 0;
}

int L_KdNodo_cmpDimension(const void *kdNodo_1, const void *kdNodo_2)
{
	double v1 = ((L_KdNodo *)kdNodo_1)->descr->c.vector[((L_KdNodo *)kdNodo_1)->dim];
	double v2 = ((L_KdNodo *)kdNodo_2)->descr->c.vector[((L_KdNodo *)kdNodo_2)->dim];
	return (v1 > v2) - (v1 < v2);
}

struct L_KdNodo_stack_data // Argumentos "virtuales" de la funcion para crear kdtrees
{
	L_KdNodo *inic;
	L_KdNodo *ret;
	int nmedio;
	int n;
	int dir; //Numero de llamada: 0=izq, 1=der; 2=salir
};

// Esta funcion se deja aca para evitar que ensucie la class_name en el .h
int L_KdTree_calcDirMayorVarianza(L_KdNodo *inic, int n)
{
	int i, j, imax;
	int nt = inic->descr->c.nt;
	std::vector<double> sx(nt), sxx(nt);
	for (i=0; i<nt; i++)
	{
		sx[i] = 0;
		sxx[i] = 0;
		for (j=0; j<n; j++)
		{
			sx[i] += inic[j].descr->c.vector[i];
			sxx[i] += inic[j].descr->c.vector[i]*inic[j].descr->c.vector[i];
		}
	}
	imax = 0;
	for (i=0; i<nt; i++)
	{
		sxx[i] = sxx[i] / n  - (sx[i]/n * sx[i]/n);
		if (sxx[i] > sxx[imax])
			imax = i;
	}
	return imax;
	
}

L_KdNodo *L_KdNodo::creaNoRec(L_KdNodo *arr, int nTot)
{
	int nit = 0;
	int i, dim;
	int nmedio;
	L_KdNodo_stack_data st[1000]; // Para eliminar recursividad (hasta donde se puede)
	int ist = 0;
	int nagr = 0;
	if (nTot==1)
	{
		L_KdNodo *ret = new L_KdNodo;
		*ret = arr[0];
		ret->dim = -1; // Es hoja
		return ret;
	}
	st[ist].inic = arr;
	st[ist].n =nTot;
	st[ist].ret = NULL;
	st[ist].nmedio = 0;
	st[ist].dir = 0;

	while (ist >= 0)
	{
		nit++;
		if (ist >= 900)
			printf("Problemas\n");
		if (st[ist].dir == 0)
		{
			if (st[ist].n == 1) // Nodo hoja
			{
				st[ist].ret = new L_KdNodo;
				*st[ist].ret = st[ist].inic[0];
				nagr++;
				st[ist].ret->dim = -1; // Es hoja
				st[ist].dir = 3;
				ist--; // "return"
				continue;
			}

			// Nodo interno
			//
			// Ordenar el arreglo segun la direccion de mayor variance
			dim = L_KdTree_calcDirMayorVarianza(st[ist].inic, st[ist].n);
			throw_L_ArgException_if(dim==-1, "L_KdTree_creaNoRec() : nodos repetidos");
			for (i=0; i<st[ist].n; i++)
				st[ist].inic[i].dim = dim;
			qsort(st[ist].inic, st[ist].n, sizeof(L_KdNodo), &L_KdNodo_cmpDimension); // Esto toma O(n*log(n)), se hace log(n) veces
			//
			// Tratar de buscar un punto de division, la idea es que inic[nmedio-1] != inic[nmedio]
			// Primero hacia la izquierda
			nmedio = st[ist].n/2;
			while (nmedio > 0 && st[ist].inic[nmedio-1].descr->c.vector[dim] == st[ist].inic[nmedio].descr->c.vector[dim])
				nmedio--;
			// Luego hacia la derecha
			if (nmedio == 0)
			{
				nmedio = st[ist].n/2;
				while (nmedio < st[ist].n && st[ist].inic[nmedio-1].descr->c.vector[dim] == st[ist].inic[nmedio].descr->c.vector[dim])
					nmedio++;
			}
			throw_L_ArgException_if(nmedio <=0 || nmedio >= st[ist].n || st[ist].inic[nmedio-1].descr->c.vector[dim] == st[ist].inic[nmedio].descr->c.vector[dim], "L_KdTree_creaNoRec() : vector repetido");
			//
			// Crear nodo de division, dividir inic[0.__n-1] -> inic[0...nmedio-1] y inic[nmedio.__n-1]
			L_KdNodo *ret = new L_KdNodo;
			*ret = st[ist].inic[nmedio];
			ret->descr = NULL;
			ret->dim = dim; // Es interior
			ret->umbral = st[ist].inic[nmedio-1].descr->c.vector[dim];

			// Empujar los datos al stack
			st[ist].ret = ret;
			st[ist].nmedio = nmedio;
			st[ist].dir = 1;

			st[ist+1].inic = st[ist].inic;
			st[ist+1].n = nmedio;
			st[ist+1].ret = NULL;
			st[ist+1].nmedio = 0;
			st[ist+1].dir = 0;
			ist++; // "call"
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
			ist++; // "call"
			continue;
		}
		else if (st[ist].dir == 2)
		{
			st[ist].ret->der = st[ist+1].ret;
			st[ist].dir = 3;
			ist--; // "return"
			continue;
		}
		throw_L_ArgException_if(st[ist].dir > 2, "L_KdTree_creaNoRec() : simulacion de recursion incorrecta");
	}
	throw_L_ArgException_if(nagr != nTot, "L_KdNodo::creaNoRec() : descriptores perdidos");
	return st[0].ret; // Decia inic, que mongo...
}


void L_KdNodo::destruyeNoRec()
{
	int si[100];
	L_KdNodo *sp[100];
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
			n++; // "call delete(izq)"
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
			n++; // "call delete(der)"
			break;
		case 2:
			si[n] = 3;
			if (sp[n]->der != NULL)
				delete sp[n]->der;
			if (sp[n]->descr != NULL)
				delete sp[n]->descr;
			n--; // "return"
		}
	} while (n>=0);
}

int L_KdNodo::cmpRepetidos(const void *d1, const void *d2)
{
	const L_KdNodo *kdes1=(const L_KdNodo *)d1, *kdes2=(const L_KdNodo *)d2;
	const L_Descriptor *des1 = &kdes1->descr->c, *des2 = &kdes2->descr->c;
	int i;
	for (i=0; i<des1->nt; i++)
	{
		if (des1->vector[i] > des2->vector[i])
			return 1;
		else if (des1->vector[i] < des2->vector[i])
			return -1;
	}
	return 0;
}

bool L_KdTreeBBF::verifyCoherenceOfObject(L_KdNodo *nodoRaiz, L_KdNodo *nodoRec)
{
	L_KdNodo *ptr;
	if (nodoRaiz==NULL)
		nodoRaiz=root;
	if (nodoRec==NULL)
		nodoRec=root;
	if (nodoRec->descr!=NULL)
	{
		if (nodoRec->izq!=NULL || nodoRec->der!=NULL)
			return false;
		// buscarme a mi mismo
		for (ptr=nodoRaiz; ptr->descr==NULL;)
		{
			if (nodoRec->descr->c.vector[ptr->dim] > ptr->umbral)
				ptr=ptr->der;
			else
				ptr=ptr->izq;
		}
		if (ptr!=nodoRec)
			return false;
		if (ptr->descr->c.nt < 0 || ptr->descr->c.nt > 10000)
			return false;
	}
	else
	{
		if (nodoRec->izq==NULL || nodoRec->der==NULL)
			return false;
		if (!verifyCoherenceOfObject(nodoRaiz,nodoRec->izq))
			return false;
		if (!verifyCoherenceOfObject(nodoRaiz,nodoRec->der))
			return false;
	}
	return true;
}

L_KdTreeBBF::L_KdTreeBBF() : paramKDTBBF("KDTBBF",0)
{
	root=NULL;
	nHojas=0;
	nElim=0;
	usePercentage=true; // Indica cual de las proximas dos opciones se usara para limitar la busqueda
	nCompMax=200; // cantidad maxima de comparaciones
	porcRev=0.25; // % de la base de datos que se revisa
	elimRepeatedElements = true; // Eliminar los nodos repetidos
	paramKDTBBF.addFrom("usePercentage",&usePercentage);
	paramKDTBBF.addFrom("nCompMax",&nCompMax);
	paramKDTBBF.addFrom("porcRev",&porcRev);
	paramKDTBBF.addFrom("elimRepeatedElements",&elimRepeatedElements);
}

bool L_KdTreeBBF::crea_noRec(L_DescriptorLista &desL)
{
	L_Array<L_KdNodo> arr((int)desL.size());

	L_DescriptorNodo *ptr;
	int i;

	if (root!=NULL)
	{
		printf("Error en L_KdTreeBBF::createFrom(): el object ya existe\n");
		return false;
	}

	if (desL.size()==0)
		return false;
	if (desL.size()==1)
	{
		root=new L_KdNodo();
		root->descr=desL.root;
		desL.root = NULL;  // Error: faltaba sacar el descriptor de desL
		desL.pult = &desL.root;
		desL.n = 0;
		return true;
	}

#ifdef L_BASIC_DEBUG
	if (desL.size()<=0)
		L_hard_shutdown("L_KdTreeBBF::createFrom : non-positive size allocation");
#endif
	// Mover los L_Node<L_Descriptor> de la lista al arreglo
	ptr=desL.root;
	for (i=0; i<(int)desL.size(); i++)
	{
		arr[i].descr=ptr; // Copia el object
		ptr=ptr->sig;
	}
	desL.root = NULL;
	desL.pult = &desL.root;
	desL.n = 0;

	nElim = arr.sortIsolatingRepeatedAtEnd(&L_KdNodo::cmpRepetidos); // Error: decia L_KdNodo3D::cmpRepetidos
	nHojas = arr.size() - nElim;

	// Recolectar nodos no absorbidos
	for (i=arr.size() - nElim; i<(int)arr.size(); i++)
	{
		*desL.pult = arr[i].descr;
		desL.pult = &(*desL.pult)->sig;
		*desL.pult = NULL;
		desL.n++;
		arr[i].descr = NULL;
	}

	arr.resize(nHojas);

	if (elimRepeatedElements)
		desL.clear();

	// Crear kdtree
	if (arr.size() == 0)
		return false;
	root = L_KdNodo::creaNoRec(arr.data(), arr.size());

	return true;
}


bool L_KdTreeBBF::crea_rec(L_DescriptorLista &desL)
{
	L_Array<L_DescriptorRefOrd> arr;

	L_DescriptorNodo *ptr;
	int i;

	if (root!=NULL)
	{
		printf("Error en L_KdTreeBBF::createFrom(): el object ya existe\n");
		return false;
	}

	if (desL.size()==0)
		return false;
	if (desL.size()==1)
	{
		root=new L_KdNodo();
		root->descr=desL.root;
		return true;
	}

#ifdef L_BASIC_DEBUG
	if (desL.size()<=0)
		L_hard_shutdown("L_KdTreeBBF::createFrom : non-positive size allocation");
#endif
	arr.resize((int)desL.size());
	ptr=desL.root;
	for (i=0; i<(int)desL.size(); i++)
	{
		arr[i].nodo=ptr;
		arr[i].asign = false;
		ptr=ptr->sig;
	}
	// NUEVO
	nElim = arr.sortIsolatingRepeatedAtEnd(&L_DescriptorRefOrd::cmpVect); // Mandar punteros con descriptores repetidos al final, dan problemas

	root=new L_KdNodo();
	// Absorber nodos de desL en forma ordenada
	if (root->creaRec(arr.data(), arr.size()-nElim, this)==false)
	{
		destroy();
		throw L_ArgException();
		return false;
	}

	//if (verifyCoherenceOfObject() == false)
	//	printf("T_T");

	desL.root = NULL;
	desL.pult = &desL.root;
	desL.n=0;
	// Recolectar nodos no absorbidos
	for (i=0; i<(int)arr.size(); i++)
	{
		if (arr[i].asign==false)
		{
			nElim++; // Nodo no fue absorbido por el tree
			*desL.pult = arr[i].nodo;
			desL.pult = &(*desL.pult)->sig;
			*desL.pult = NULL;
			desL.n++;
		}
		else
			nHojas++; // Nodo fue absorbido por el tree
	}
	return true;
}


bool L_KdTreeBBF::createFrom(L_DescriptorLista &desL)
{
	//return crea_rec(desL);
	return crea_noRec(desL);
}


// addFrom() no soporta la opcion delete_or_isolate_repeated=true. Hace lo mismo que con delete_or_isolate_repeated=false.
bool L_KdTreeBBF::addFrom(L_DescriptorLista &desL) // addFrom un nodo a un tree ya creado => esto desbalancea el tree
{
	L_KdNodo *ptrN;
	L_DescriptorNodo *ptrD;

	if (root==NULL) // tree en blanco
	{
		return createFrom(desL); // Crear tree balanceado del inicio
	}
	ptrD=desL.root;

	// tree absorbe nodos de desL en forma ordenada
	while(ptrD!=NULL)
	{
		ptrN=root;
		// Agregar el nodo ptrD al tree
		while (ptrN!=NULL)
		{
			if (ptrN->descr!=NULL) // Se llega a un nodo hoja. Hay que transformarlo en nodo de division
			{
				// Buscar dimension de maxima variance
				double var, varmax=0;
				int i, imax=0;
				for (i=0; i<ptrD->c.nt; i++)
				{
					var=ptrD->c.vector[i]-ptrN->descr->c.vector[i];
					var*=var;
					if (var > varmax)
					{
						varmax=var;
						imax=i;
					}
				}
				if (varmax==0) // Este descriptor ya existe en el tree. Eliminarlo.
				{
					L_DescriptorNodo *tempD=ptrD;
					ptrD=ptrD->sig;
					tempD->sig=NULL;
					delete tempD;
				}
				// Transformar nodo hoja en nodo de division con dos children nodo hoja
				ptrN->dim=imax;
				if (ptrD->c.vector[imax]<ptrN->descr->c.vector[imax])
				{
					ptrN->umbral=ptrD->c.vector[imax];
					ptrN->izq=new L_KdNodo();
					ptrN->izq->descr=ptrD;
					ptrN->der=new L_KdNodo();
					ptrN->der->descr=ptrN->descr;
					ptrN->descr=NULL;
				}
				else
				{
					ptrN->umbral=ptrN->descr->c.vector[imax];
					ptrN->izq=new L_KdNodo();
					ptrN->izq->descr=ptrN->descr;
					ptrN->der=new L_KdNodo();
					ptrN->der->descr=ptrD;
					ptrN->descr=NULL;
				}
			}
			else // Se llega a un nodo de division
			{
				if (ptrD->c.vector[ptrN->dim]<=ptrN->umbral) // Ir hacia la izquierda
					ptrN=ptrN->izq;
				else // Ir hacia la derecha
					ptrN=ptrN->der;
			}
		}
	}
	// Eliminar la referencia de desL a los nodos absorbidos
	desL.root=NULL;
	desL.pult=&desL.root;
	desL.n=0;

	return true;
}

bool L_KdTreeBBF::optimize()
{
	throw L_NonImplementedException();
	return false;
}

bool L_KdTreeBBF::destroy()
{
	if (root!=NULL)
	{
		//root->destroyRec();
		root->destruyeNoRec();
		delete root;
		root=NULL;
	}
	nHojas=0;
	nElim=0;
	return true;
}

L_ParamManagerLocal *L_KdTreeBBF::pideParams()
{
	return &paramKDTBBF;
}

//! Devuelve los n nodos del árbol que son más parecidos a orig en {dest[0],dest[1],...}.
void L_KdTreeBBF::findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	int nComp=0;
	L_BBFPriorLista lisPrio;
	L_KdNodo *ptr;
	L_Descriptor *des;
	double dif;
	double difTemp;
	L_Descriptor *desTemp;
	int i;
	int nBus;

	if (usePercentage)
		nBus=(int)(1+floor(porcRev*nHojas));
	else
		nBus=nCompMax;

	if (nBus>nHojas)
		nBus=nHojas;

	for (i=0; i<n; i++)
	{
		dest[i]=NULL;
		dist[i]=L_BBFPRIORNODO_DISTINF; // número grande
	}

	if (root==NULL)
		return;
	else if (root->descr!=NULL)
	{
		dist[0]=orig.distEuclidCuadr_mismoTipo(root->descr->c);
		dest[0]=&root->descr->c;
		return;
	}

	dif=orig.vector[root->dim]-root->umbral;

	if (dif>0)
	{
		lisPrio.addFrom(*root, *root->der, orig, 0); // ir primero por la derecha, distancia puesta a 0
		lisPrio.addFrom(*root, *root->izq, orig);
	}
	else
	{
		lisPrio.addFrom(*root, *root->izq, orig, 0); // ir primero por la izquierda, distancia puesta a 0
		lisPrio.addFrom(*root, *root->der, orig);
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
			if (ptr->descr!=NULL) // ptr es nodo hoja
				break;
			else // ptr es nodo de division
			{
				dif=orig.vector[ptr->dim]-ptr->umbral;
				if (dif>0) // bajar por la derecha
				{
					lisPrio.addFrom(*ptr, *ptr->izq, orig);
					ptr=ptr->der;
				}
				else // bajar por la izquierda
				{
					lisPrio.addFrom(*ptr, *ptr->der, orig);
					ptr=ptr->izq;
				}
			}
		}
		des=&ptr->descr->c;
		dif=sqrt(des->distEuclidCuadr_mismoTipo(orig));
		nComp++;
		for (i=0; i<n; i++)
		{
			if (dest[i]==NULL)
			{
				dest[i]=des;
				dist[i]=dif;
				break;
			}
			if (dif<dist[i])
			{
				difTemp=dif;
				desTemp=des;
				dif=dist[i];
				des=dest[i];
				dist[i]=difTemp;
				dest[i]=desTemp;
			}
		}
	}
}

//! Devuelve los n nodos del árbol que son más parecidos a orig en {dest[0],dest[1],...}. Cada dest[i] debe pertenecer a un object único.
void L_KdTreeBBF::findNearest_otherObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	int nComp=0;
	L_BBFPriorLista lisPrio;
	L_KdNodo *ptr;
	L_Descriptor *des;
	double dif;
	int i;
	bool objYaEsta;
	int nBus;

	if (usePercentage)
		nBus=(int)(1+floor(porcRev*nHojas));
	else
		nBus=nCompMax;
	if (nBus>nHojas)
		nBus=nHojas;

	for (i=0; i<n; i++)
	{
		dist[i]=L_BBFPRIORNODO_DISTINF; // número grande
		dest[i]=NULL;
	}

	if (root==NULL)
		return;
	else if (root->descr!=NULL)
	{
		dist[0]=orig.distEuclidCuadr_mismoTipo(root->descr->c);
		dest[0]=&root->descr->c;
	}

	dif=orig.vector[root->dim]-root->umbral;

	if (dif>0)
	{
		lisPrio.addFrom(*root, *root->der, orig); // ir primero por la derecha
		lisPrio.addFrom(*root, *root->izq, orig);
	}
	else
	{
		lisPrio.addFrom(*root, *root->izq, orig); // ir primero por la izquierda
		lisPrio.addFrom(*root, *root->der, orig);
	}

	while (true)
	{
		if (nComp>=nBus) // se acabó la búsqueda
			break;
		ptr=lisPrio.popNodoMasCercano();
		if (ptr==NULL) // se acabó la búsqueda
			break;
		// Bajar por el árbol hasta llegar a una hoja => des
		while (true)
		{
			if (ptr->descr!=NULL) // ptr es nodo hoja
			{
				des=&ptr->descr->c;
				break;
			}
			else // ptr es nodo de division
			{
				dif=orig.vector[ptr->dim]-ptr->umbral;
				if (dif>0) // bajar por la derecha
				{
					lisPrio.addFrom(*ptr, *ptr->izq, orig);
					ptr=ptr->der;
				}
				else // bajar por la izquierda
				{
					lisPrio.addFrom(*ptr, *ptr->der, orig);
					ptr=ptr->izq;
				}
			}
		}
		dif=sqrt(des->distEuclidCuadr_mismoTipo(orig));
		nComp++;
		objYaEsta=false;
		for (i=0; i<n; i++)
		{
			if (dest[i]!=NULL && dest[i]->nobj==des->nobj)
			{
				objYaEsta=true;
				break;
			}
		}
		if (objYaEsta)
		{
			int iOtro=i;
			if (dif<dist[iOtro])
			{
				dist[iOtro]=L_BBFPRIORNODO_DISTINF; // Este muere seguro.
				if (dif<dist[n-1])
				{
					double difTemp;
					L_Descriptor *desTemp;
					for (i=0; i<=iOtro; i++)
					{
						if (dif<dist[i])
						{
							difTemp=dif;
							desTemp=des;
							dif=dist[i];
							des=dest[i];
							dist[i]=difTemp;
							dest[i]=desTemp;
						}
					}
				}
			}
		}
		else
		{
			if (dif<dist[n-1])
			{
				double difTemp;
				L_Descriptor *desTemp;
				for (i=0; i<n; i++)
				{
					if (dif<dist[i])
					{
						difTemp=dif;
						desTemp=des;
						dif=dist[i];
						des=dest[i];
						dist[i]=difTemp;
						dest[i]=desTemp;
					}
				}
			}
		}
	}
}

// addFrom los nodos en primer lugar por distancia al objetivo, y en segundo lugar por tiempo de llegada
void L_BBFPriorLista::addFrom(L_KdNodo &nododiv, L_KdNodo &nododir, const L_Descriptor &orig, double dist) // addFrom nodo en orden de mas cercania a orig
{
	L_BBFPriorNodo *nuevo;
	L_BBFPriorNodo **pptr;

	if (dist<-0.5)
		L_SETABS(dist, orig.vector[nododiv.dim]-nododiv.umbral);
	if (dist>cotaDist)
		return;

	new_L_BBFPriorNodo(nuevo,(*this));
	nuevo->dist=dist;
	nuevo->dir=&nododir;
	pptr=&root;
	while (*pptr!=NULL && (*pptr)->dist<=dist)
		pptr=&(*pptr)->sig;
	nuevo->sig=*pptr;
	*pptr=nuevo;
	n++;
}

L_KdNodo *L_BBFPriorLista::popNodoMasCercano() // devuelve nodo sin memoria propia
{
	L_BBFPriorNodo *pri;
	L_KdNodo *ret;

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
	delete_L_BBFPriorNodo(pri, (*this));
	n--;
	return ret;
}

bool L_BusqSecDescr::createFrom(L_DescriptorLista &desL)
{
	L_DescriptorNodo *des;
	L_DescriptorPtr desPtr;
	if (lista.size()!=0)
	{
		printf("Error en L_BusqSecDescr::createFrom(): el object ya existe\n");
		return false;
	}
	for (des=desL.root; des!=NULL; des=des->sig)
	{
		desPtr.ptr=des;
		lista.push_back(desPtr);
	}
	if (delete_or_isolate_repeated)
		lista.sort(2, &L_Descriptor::qCmpMula_nodo);
	else
		lista.sort(1, &L_Descriptor::qCmpMula_nodo);
	return true;
}

bool L_BusqSecDescr::addFrom(L_DescriptorLista &desL)
{
	L_DescriptorNodo *des;
	L_DescriptorPtr desPtr;
	for (des=desL.root; des!=NULL; des=des->sig)
	{
		desPtr.ptr=des;
		lista.push_back(desPtr);
	}
	if (delete_or_isolate_repeated)
		lista.sort(1, &L_Descriptor::qCmpMula_nodo);
	else
		lista.sort(1, &L_Descriptor::qCmpMula_nodo);
	return true;
}

bool L_BusqSecDescr::destroy()
{
	lista.clear(); // No destroy los descriptores reales
	return true;
}

bool L_BusqSecDescr::optimize()
{
	return false;
}

void L_BusqSecDescr::findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	L_DescriptorPtrNodo *ptr;
	double d;
	int i,j;
	for (i=0; i<n; i++)
		dest[i]=NULL;
	for (ptr=lista.root; ptr!=NULL; ptr=ptr->sig)
	{
		for (i=0; i<n; i++)
		{
			if (dest[i]==NULL)
			{
				dest[i]=&ptr->c.ptr->c;
				dist[i]=orig.distEuclidCuadr_mismoTipo(ptr->c.ptr->c);
			}
			else
			{
				d=orig.distEuclidCuadr_mismoTipo(ptr->c.ptr->c);
				if (d<dist[i])
				{
					for (j=n-1; j>i; j--)
					{
						dest[j]=dest[j-1];
						dist[j]=dist[j-1];
					}
					dest[i]=&ptr->c.ptr->c;
					dist[i]=d;
					break;
				}
			}
		}
	}
	return;
}

void L_BusqSecDescr::findNearest_otherObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	L_DescriptorPtrNodo *ptr;
	double d;
	int i,j;
	for (i=0; i<n; i++)
		dest[i]=NULL;
	for (ptr=lista.root; ptr!=NULL; ptr=ptr->sig)
	{
		for (i=0; i<n; i++)
		{
			if (dest[i]==NULL)
			{
				dest[i]=&ptr->c.ptr->c;
				dist[i]=orig.distEuclidCuadr_mismoTipo(ptr->c.ptr->c);
			}
			else
			{
				d=orig.distEuclidCuadr_mismoTipo(ptr->c.ptr->c);
				if (d<dist[i])
				{
					for (j=i; j<n; j++)
					{
						if (dest[j]->nobj==ptr->c.ptr->c.nobj)
							break;
					}
					for (; j>i; j--)
					{
						dest[j]=dest[j-1];
						dist[j]=dist[j-1];
					}
					dest[i]=&ptr->c.ptr->c;
					dist[i]=d;
					break;
				}
				else
				{
					if (dest[i]->nobj==ptr->c.ptr->c.nobj)
					{
						break;
					}
				}
			}
		}
	}
	return;
}

L_BusqSecDescr::boolVBC L_BusqSecDescr::verifBusqCorrecta(const L_Descriptor &orig, const L_Descriptor &prim, const L_Descriptor &seg)
{
	L_DescriptorPtrNodo *ptr;
	double min1=100000.0, min2=100000.0; // valen "infinito"
	double dist;
	for (ptr=lista.root; ptr!=NULL; ptr=ptr->sig)
	{
		dist=orig.distEuclidCuadr_mismoTipo(ptr->c.ptr->c);
		if (dist<min2)
		{
			if (dist<min1)
			{
				min2=min1;
				min1=dist;
			}
			else
			{
				min2=dist;
			}
		}
	}
	if (orig.distEuclidCuadr_mismoTipo(prim) == min1)
	{
		if (orig.distEuclidCuadr_mismoTipo(seg) == min2)
		{
			return L_BusqSecDescr::sin_error;
		}
		return L_BusqSecDescr::error_segundo;
	}
	if (orig.distEuclidCuadr_mismoTipo(seg) == min2)
	{
		return L_BusqSecDescr::error_primero;
	}
	return L_BusqSecDescr::error_ambos;
}

bool L_BusqPermPiv::createFrom(L_DescriptorLista &desL)
{
	int i, nTodos, nSorteo;
	L_DescriptorNodo *des;
	L_DescriptorPtrLista li;
	L_DescriptorPtr desPtr;
	L_DescriptorPtrNodo *nodo;
	std::vector<int> todos, sorteo;
	if (bd.size()!=0)
	{
		printf("Error en L_BusqPermPiv::createFrom(): ya existe el object\n");
		return false; // No se puede crear, ya existe
	}
	if (desL.size()==0)
		return false; // No tiene sentido. Revisar que hacer aca
	// Eliminar nodos repetidos
	for (des=desL.root; des!=NULL; des=des->sig)
	{
		desPtr.ptr=des;
		li.push_back(desPtr);
	}
	if (delete_or_isolate_repeated)
		li.sort(2, &L_Descriptor::qCmpMulaPtr_nodo);
	else
		li.sort(1, &L_Descriptor::qCmpMulaPtr_nodo);
	// Absorber lista li
	nBd=(int)li.size();
#ifdef L_BASIC_DEBUG
	if (nBd<=0)
		L_hard_shutdown("L_BusqPermPiv::createFrom : intento de crear object de tamano nBd = zero");
#endif
	bd.resize(nBd);
	for (i=0,nodo=li.root; i<nBd; i++,nodo=nodo->sig)
	{
		bd[i].des=nodo->c.ptr;
		bd[i].perm.resize(nPiv);
	}
	// Elegir pivotes. Los indices de los pivotes quedan en sorteo[]
	nTodos=nBd;
	nSorteo=0;
	todos.resize(nTodos);
#ifdef L_BASIC_DEBUG
	if (nPiv<=0)
		L_hard_shutdown("L_BusqPermPiv::createFrom : intento de crear object de tamano nPiv = zero");
#endif
	sorteo.resize(nPiv);
	for (i=0; i<nTodos; i++)
		todos[i]=i;
	while (nBd-nTodos < nPiv)
	{
		i=(int)( rand()%nTodos );
		throw_L_ArgException_if(i<0 || i>=nTodos, "L_BusqPermPiv::createFrom : Error de rango");
		sorteo[nSorteo++]=i;
		todos[i]=todos[nTodos-1];
		nTodos--;
	}
	// Generar pivotes
	piv.resize(nPiv);
	for (i=0; i<nPiv; i++)
		piv[i]=bd[sorteo[i]].des;
	// Calcular distancias de base de datos a pivotes y permutaciones
	dist.resize(nPiv);
	for (i=0; i<nBd; i++)
	{
		calcPerm(bd[i].des->c, bd[i]);
	}
	// init memoria temporal para futuro calculo de regla de Spearman
	temp.resize(nPiv);
	return true;
}

void L_BusqPermPiv::calcPerm(const L_Descriptor &des, L_PermPiv &permObj)
{
	int i;
	for (i=0; i<nPiv; i++)
	{
		dist[i].ind=i;
		dist[i].d=piv[i]->c.distEuclidCuadr_mismoTipo(des);
	}
	qsort(&(dist[0]), nPiv, sizeof(dist_str), &cmpDist_str);
	for (i=0; i<nPiv; i++)
		permObj.perm[i]=dist[i].ind;
}

int L_BusqPermPiv::cmpDist_str(const void *a, const void *b)
{
	return ( ((dist_str *)a)->d > ((dist_str *)b)->d )
		- ( ((dist_str *)a)->d < ((dist_str *)b)->d );
}

int L_BusqPermPiv::calcSpearman(const L_PermPiv &p1, const L_PermPiv &p2)
{
	int i, total=0;
	for(i=0; i<nPiv; i++)
		temp[p1.perm[i]] = i;   // compone la primera secuencia a su forma canonica
	for(i=0; i<nPiv; i++)
		total += abs((i-temp[p2.perm[i]]));
	return total;
}

int L_BusqPermPiv::cmpPermDist(const void *a, const void *b)
{
	return ( ((L_PermPiv *)a)->dist > ((L_PermPiv *)b)->dist )
		- ( ((L_PermPiv *)a)->dist < ((L_PermPiv *)b)->dist );
}

void L_BusqPermPiv::ordenaPiv(const L_Descriptor &orig)
{
	int i;
	L_PermPiv origPerm;
	origPerm.perm.resize(nPiv);
	calcPerm(orig, origPerm);
	for (i=0; i<nBd; i++)
		bd[i].dist=calcSpearman(bd[i],origPerm);
	qsort(&(bd[0]), nBd, sizeof(L_PermPiv), &cmpPermDist);
}

void L_BusqPermPiv::findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	int i, j;
	L_Descriptor *de;
	double di;
	int nBus;
	L_Descriptor *detemp;
    double ditemp;

	ordenaPiv(orig);
	for (i=0; i<n; i++)
		dest[i]=NULL; // Muy grande
	if (usePercentage)
		nBus=(int)(1+floor(porcRev*nBd));
	else
		nBus=nCompMax;

	if (nBus>nBd)
		nBus=nBd;

	for (i=0; i<nBus; i++)
	{
		de=&bd[i].des->c;
		di=bd[i].des->c.distEuclidCuadr_mismoTipo(orig);
		for (j=0; j<n; j++)
		{
			if (dest[j]==NULL)
			{
				dest[j]=de;
				dist[j]=di;
				break;
			}
			if (di<dist[j])
			{
				detemp=dest[j];
				ditemp=dist[j];
				dest[j]=de;
				dist[j]=di;
				de=detemp;
				di=ditemp;
			}
		}
	}
}

void L_BusqPermPiv::findNearest_otherObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	int i, j;
	L_Descriptor *de;
	double di;
	int nBus;
	L_Descriptor *detemp;
    double ditemp;

	ordenaPiv(orig);
	for (i=0; i<n; i++)
		dest[i]=NULL; // Muy grande
	if (usePercentage)
		nBus=(int)(1+floor(porcRev*nBd));
	else
		nBus=nCompMax;

	if (nBus>nBd)
		nBus=nBd;

	for (i=0; i<nBus; i++)
	{
		de=&bd[i].des->c;
		di=bd[i].des->c.distEuclidCuadr_mismoTipo(orig);
		for (j=0; j<n; j++)
		{
			if (dest[j]==NULL)
			{
				dest[j]=de;
				dist[j]=di;
				break;
			}
			if (di<dist[j])
			{
				detemp=dest[j];
				ditemp=dist[j];
				dest[j]=de;
				dist[j]=di;
				de=detemp;
				di=ditemp;
			}
			if (de->nobj==dest[j]->nobj)
				break;
		}
	}
}

bool L_BusqPermPiv::optimize()
{
	return false;
}

L_ParamManagerLocal *L_BusqPermPiv::pideParams()
{
	return &paramBusPP;
}

#ifdef FLANN_H
void L_Flann::fijaParamB()
{
	flannParams.log_level = LOG_WARN; //LOG_INFO;
	flannParams.log_destination = NULL;
    flannParams.algorithm = KDTREE;
    flannParams.checks = 32;
    flannParams.trees = 8;
    flannParams.branching = 32;
    flannParams.iterations = 7;
    flannParams.target_precision = -1;
}
#endif

#ifdef FLANN_H
bool L_Flann::createFrom(L_DescriptorLista &desc)
{
	float speedup;
	int i, j, nt;
	if (desc.root == NULL)
		return false;
	nt = desc.root->c.nt;
	desc.copia_a_arreglo(desArr);
	desc.clear();
	data.resize(desArr.size() * nt);
	data.adjustMemDown();
	for (i=0; i<desArr.size(); i++)  // = desArr.size()
		for (j=0; j<nt; j++)
			data[i*nt+j] = (float)desArr[i].vector[j];
	index = flann_build_index(data.elem, desArr.size(), nt, &speedup, &flannParams);
	if (index == NULL)
		return false;
	return true;
}
#endif

#ifdef FLANN_H
bool L_Flann::destroy()
{
	if (index != NULL)
	{
		flann_free_index(index, &flannParams);
		index = NULL;
	}
	desArr.clear();
	data.clear();
	return true;
}
#endif

#ifdef FLANN_H
void L_Flann::findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	std::vector<float> testset(orig.nt);
	std::vector<int> indices(n);
	std::vector<float> dists(n);
	int i, checks;
	for (i=0; i<orig.nt; i++)
		testset[i] = (float)orig.vector[i];
	for (i=0; i<n; i++)
		dest[i] = NULL;
	if (index==NULL)
		return;

	checks = flannParams.checks; // 50 o flannParams.checks
	flann_find_nearest_neighbors_index(index, testset.elem, 1, indices.elem, dists.elem, n, checks, &flannParams);
	for (i=0; i<n; i++)
	{
		dest[i] = &desArr[indices[i]];
		dist[i] = dists[i];
	}
}
#endif

#ifdef FLANN_H
void L_Flann::findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	findNearest_anyObject(orig, n, dest, dist); // No hay implementacion especial
}
#endif

#ifdef FLANN_H
void L_Flann::encuentraMasCercanos(std::vector<L_Descriptor> &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
{
	int i, j, checks;
	for (i=0; i<n*orig.size(); i++)
		dest[i] = NULL;
	if (index==NULL || orig.elem == NULL)
		return;
	int nt = orig[0].nt;
	std::vector<float> testset(orig.size()*nt);
	std::vector<int> indices(orig.size()*n);
	std::vector<float> dists(orig.size()*n);
	for (i=0; i<orig.size(); i++)
		for (j=0; j<nt; j++)
			testset[i*nt+j] = (float)orig[i].vector[j];

	checks = flannParams.checks; // 50 o flannParams.checks
	flann_find_nearest_neighbors_index(index, testset.elem, orig.size(), indices.elem, dists.elem, n, flannParams.checks, &flannParams);

	for (i=0; i<orig.size(); i++)
	{
		for (j=0; j<n; j++)
		{
			dest[i*n+j] = &desArr[indices[i*n+j]];
			dist[i*n+j] = dists[i*n+j];
		}
	}
}
#endif


void L_CalceLista::genLineas(L_ShapeArray &lins, int nImRef)
{
	L_CalceNodo *ptr;
	int nir, nip = 0; // A la imagen de test se le asign el numero global de imagen 0

	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		nir = ptr->c.dRef->nobj+1; // A la imagen de referencia nobj se le asign el numero global de imagen nobj+1
		if (nImRef == -1 || nImRef == nir)
			lins._drawShape(nir, nip, (int)ptr->c.dRef->x0, (int)ptr->c.dRef->y0, (int)ptr->c.dPru->x0, (int)ptr->c.dPru->y0);
	}
}

double L_CalceLista::calcAngProm()
{
	double x=0, y=0;
	L_CalceNodo *ptr;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		x+=cos(ptr->c.dRef->ang - ptr->c.dPru->ang);
		y+=sin(ptr->c.dRef->ang - ptr->c.dPru->ang);
	}
	if (y==0 && x==0)
		return 0;
	else
		return atan2(y,x);
}

bool L_CalceLista::revisaCoherencia()
{
	double sx0=0, sy0=0, sv0=0;
	L_CalceNodo *cal;
	for (cal=root; cal!=NULL; cal=cal->sig)
	{
		if (cal->c.dPru == NULL || cal->c.dRef == NULL)
			return false;
		sx0 += cal->c.dPru->x0 * cal->c.dPru->x0;
		sy0 += cal->c.dPru->y0 * cal->c.dPru->y0;
		sv0 += cal->c.dPru->vector[0] * cal->c.dPru->vector[0];
		sx0 += cal->c.dRef->x0 * cal->c.dRef->x0;
		sy0 += cal->c.dRef->y0 * cal->c.dRef->y0;
		sv0 += cal->c.dRef->vector[0] * cal->c.dRef->vector[0];
	}
	return (sx0+sy0+sv0) != 0;
}

bool L_CalceLista::restrSemiLocales(const L_Calce &central, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng)
{
	L_Array<L_Calce *> vecinosRef(nVecinos);
	L_Array<L_Calce *> vecinosPru(nVecinos);
	std::vector<double> distVecinosRef(nVecinos);
	std::vector<double> distVecinosPru(nVecinos);
	L_Array<L_Calce *> corresp(nVecinos);
	int i, ii, j, nCorresp=0;
	double dist;
	double d1, d2, d3, d4, dx1, dy1, dx2, dy2, ang1, ang2;
	L_CalceNodo *cal;

	if (nVecinos >= this->n)
		return true;
	for (i=0; i<nVecinos; i++)
	{
		distVecinosRef[i]=1e30;
		distVecinosPru[i]=1e30;
		vecinosRef[i]=NULL;
		vecinosPru[i]=NULL;
	}
	for (cal=root ; cal!=NULL ; cal=cal->sig)
	{
		if (&cal->c==&central)
			continue;
		dist=cal->c.dRef->distEuclidCuadr(*central.dRef);
		if (dist < distVecinosRef[nVecinos-1])
		{
			for (i=nVecinos-2; i>=0; i--)
			{
				if (dist > distVecinosRef[i])
				{
					vecinosRef[i+1]=&cal->c;
					distVecinosRef[i+1]=dist;
					break;
				}
				else if (i==0)
				{
					vecinosRef[i]=&cal->c;
					distVecinosRef[i]=dist;
					break;
				}
				else
				{
					vecinosRef[i+1]=vecinosRef[i];
					distVecinosRef[i+1]=distVecinosRef[i];
				}
			}
		}
		dist=cal->c.dPru->distEuclidCuadr(*central.dPru);
		if (dist < distVecinosPru[nVecinos-1])
		{
			for (i=nVecinos-2; i>=0; i--)
			{
				if (dist > distVecinosPru[i])
				{
					vecinosPru[i+1]=&cal->c;
					distVecinosPru[i+1]=dist;
					break;
				}
				else if (i==0)
				{
					vecinosPru[i]=&cal->c;
					distVecinosPru[i]=dist;
					break;
				}
				else
				{
					vecinosPru[i+1]=vecinosPru[i];
					distVecinosPru[i+1]=distVecinosPru[i];
				}
			}
		}
	}
	// Encontrar los calces que tienen ambos puntos de interés cerca del central, los que aqui son llamados correspondencias
	for (i=0; i<nVecinos; i++)
	{
		if (vecinosRef[i]==NULL) // No deberia pasar, pero por si aca...
			continue;
		for (j=0; j<nVecinos; j++)
		{
			if (vecinosPru[j]==NULL) // No deberia pasar, pero por si aca...
				continue;
			if (vecinosRef[i]==vecinosPru[j])
			{
				corresp[nCorresp]=vecinosRef[i];
				nCorresp++;
				break;
			}
		}
	}
	if (nCorresp<minCorresp)
		return false; // Muy pocas correspondencias
	// Verificar que distancias y angulos entre correspondencias sean parecidos
	if (verifEscalaAngs)
	{
		for (i=0; i<nCorresp; i++)
		{
			// Se revisan los triangulos {vecinos[i], central, vecinos[i+1]} en la imagen de referencia y de test
			// triangulo en imagen de referencia
			ii = (i+1) % nCorresp;
			dx1=corresp[i]->dRef->x0-central.dRef->x0;
			dy1=corresp[i]->dRef->y0-central.dRef->y0;
			dx2=corresp[ii]->dRef->x0-central.dRef->x0;
			dy2=corresp[ii]->dRef->y0-central.dRef->y0;
			d1=dx1*dx1+dy1*dy1;
			d2=dx2*dx2+dy2*dy2;
			ang1=asin( L_CrossProduct2D(dx1, dy1, dx2, dy2, double) );
			// triangulo en imagen de test
			dx1=corresp[i]->dPru->x0-central.dPru->x0;
			dy1=corresp[i]->dPru->y0-central.dPru->y0;
			dx2=corresp[ii]->dPru->x0-central.dPru->x0;
			dy2=corresp[ii]->dPru->y0-central.dPru->y0;
			d3=dx1*dx1+dy1*dy1;
			d4=dx2*dx2+dy2*dy2;
			ang2=asin( L_CrossProduct2D(dx1, dy1, dx2, dy2, double) );
			// restricciones semi-locales
			if (d1==0 || d2==0 || d3==0 || d4==0)
				return false;
			if ( (d1/d2) / (d3/d4) > varEscala*varEscala || (d3/d4) / (d1/d2) > varEscala*varEscala )
				return false;
			if (fabs(ang1-ang2) > varAng)
				return false;
		}
	}
	return true;
}

int L_CalceLista::aplicaRestrSemiLocales(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng)
{
	L_CalceNodo *tmp=NULL;
	int nFallados=0;
	n=0;
	for (pult=&root; *pult!=NULL;)
	{
		if (restrSemiLocales((**pult).c, nVecinos, minCorresp, verifEscalaAngs, varEscala, varAng))
		{
			pult=&(*pult)->sig;
			n++;
		}
		else
		{
			nFallados++;
			if (corrigeLista)
			{
				tmp=*pult;
				*pult=(*pult)->sig;
				delete tmp;
			}
			else
			{
				n++;
				pult=&(*pult)->sig;
			}
		}
	}
	//printf("Fallados: %d de %d\n", nFallados, n+nFallados);
	return nFallados;
}

bool L_CalceLista::cuentaCerros(L_Calce &c1, L_Calce &c2, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin)
{
	// Idea: guardar pixeles de cada linea de sonar 1 en arreglo unidimensional
	// En cada arreglo unidimensional pasar pares de filtros IR pasabandas (para simetria) y luego contar cerros
	L_SignalDouble sPru, sRef;
	L_SignalDouble sMinPru, sMinRef, sMaxPru, sMaxRef;
	int semiancho;
	int nMaxPru=0, nMinPru=0, nMaxRef=0, nMinRef=0;
	int num1, num2;
	sRef.scanLineFromImage(*c1.dRef->imorig,c1.dRef->x0,c1.dRef->y0,c2.dRef->x0,c2.dRef->y0);
	sPru.scanLineFromImage(*c1.dPru->imorig,c1.dPru->x0,c1.dPru->y0,c2.dPru->x0,c2.dPru->y0);
	sRef.normalizeHistogram();
	sPru.normalizeHistogram();
	L_VM(sRef.s).sum_to_each_element(-L_VMC(sRef.s).mean());
	L_VM(sPru.s).sum_to_each_element(-L_VMC(sPru.s).mean());
	sRef[0]=0;
	sRef[sRef.size()-1]=0;
	sPru[0]=0;
	sPru[sPru.size()-1]=0;
	sRef.filter_IR_right(0.5);
	sPru.filter_IR_left(0.5);

	semiancho=(int)(minDistCerros/2);
	sRef.countMaxMin(nMaxRef, nMinRef, -0.1, 0.1, semiancho);
	sPru.countMaxMin(nMaxPru, nMinPru, -0.1, 0.1, semiancho);

	/*
	{
		L_ImageRGBUchar ima, imb;
		sRef.genGraphic(ima);
		sPru.genGraphic(imb);
		ima.writeImageATL("graf_ref.tif");
		imb.writeImageATL("graf_pru.tif");
	}
	*/
	num1=L_MAX(nMaxPru,nMaxRef);
	num2=L_MIN(nMaxPru,nMaxRef);

	if (num2 < minCerrosUtil)
		return false;
	if (num2 /(double) num1 < parecidoCerrosMin)
		return false;

	num1=L_MAX(nMinPru,nMinRef);
	num2=L_MIN(nMinPru,nMinRef);

	if (num2 < minCerrosUtil)
		return false;
	if (num2 /(double) num1 < parecidoCerrosMin)
		return false;

	return true;
}

double L_CalceLista::aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax)
{
	//L_CalceNodo *tmp=NULL;
	long nExitosos=0;
	long nRevisados=0;
	L_Array<L_Calce *> arr;
	int u, i, j;
	double x1, y1, x2, y2;

	arr.resize(n);
	for (i=0, pult=&root; i<n; i++, pult=&(*pult)->sig)
		arr[i]=&(*pult)->c;

	if (n*n < nLineasMax)
	{
		for (i=0; i<n-1; i++)
			for (j=i+1; j<n; j++)
			{
				x1=arr[i]->dPru->x0;
				y1=arr[i]->dPru->y0;
				x2=arr[j]->dPru->x0;
				y2=arr[j]->dPru->y0;
				if ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) > 25*25)
				{
					nExitosos+=cuentaCerros(*arr[i],*arr[j], minDistCerros, maxDistCerros, minCerrosUtil, parecidoCerrosMin)?1:0;
					nRevisados++;
				}
			}
	}
	else
	{
		for (u=0; u<nLineasMax; u++)
		{
			i=rand()%n;
			j=rand()%n;
			if (i==j)
				continue;
			x1=arr[i]->dPru->x0;
			y1=arr[i]->dPru->y0;
			x2=arr[j]->dPru->x0;
			y2=arr[j]->dPru->y0;
			if ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) > 25*25)
			{
				nExitosos+=cuentaCerros(*arr[i],*arr[j], minDistCerros, maxDistCerros, minCerrosUtil, parecidoCerrosMin)?1:0;
				nRevisados++;
			}
		}
	}
	if (nRevisados==0)
	{
		nExitosos=1; // Para evitar 0/0
		nRevisados=1;
	}
	return nExitosos/(double)nRevisados;
}

int L_CalceLista::cmpCalceXY(const void *a, const void *b)
{
	L_CalceNodo *c1, *c2;
	c1=*(L_CalceNodo **)a;
	c2=*(L_CalceNodo **)b;
	if (c1->c.dPru->x0 > c2->c.dPru->x0)
		return 1;
	if (c1->c.dPru->x0 < c2->c.dPru->x0)
		return -1;
	if (c1->c.dPru->y0 > c2->c.dPru->y0)
		return 1;
	if (c1->c.dPru->y0 < c2->c.dPru->y0)
		return -1;
	if (c1->c.dRef->x0 > c2->c.dRef->x0)
		return 1;
	if (c1->c.dRef->x0 < c2->c.dRef->x0)
		return -1;
	if (c1->c.dRef->y0 > c2->c.dRef->y0)
		return 1;
	if (c1->c.dRef->y0 < c1->c.dRef->y0)
		return -1;
	return 0;
}

int L_CalceLista::cmpCalceXY_v2(const void *a, const void *b)
{
	L_CalceNodo *c1, *c2;
	c1=*(L_CalceNodo **)a;
	c2=*(L_CalceNodo **)b;
	if (c1->c.dRef->x0 > c2->c.dRef->x0)
		return 1;
	if (c1->c.dRef->x0 < c2->c.dRef->x0)
		return -1;
	if (c1->c.dRef->y0 > c2->c.dRef->y0)
		return 1;
	if (c1->c.dRef->y0 < c1->c.dRef->y0)
		return -1;
	if (c1->c.dPru->x0 > c2->c.dPru->x0)
		return 1;
	if (c1->c.dPru->x0 < c2->c.dPru->x0)
		return -1;
	if (c1->c.dPru->y0 > c2->c.dPru->y0)
		return 1;
	if (c1->c.dPru->y0 < c2->c.dPru->y0)
		return -1;
	return 0;
}

void L_CalceLista::transformaEnMatrices(L_Matrix &uv_1, L_Matrix &uv_2)
{
	int i=0;
	uv_1.reallocate(n, 2);
	uv_2.reallocate(n, 2);
	L_CalceNodo *c;
	for (c=root; c!=NULL; c=c->sig,i++)
	{
		uv_1(i,0)=c->c.dRef->x0;
		uv_1(i,1)=c->c.dRef->y0;
		uv_2(i,0)=c->c.dPru->x0;
		uv_2(i,1)=c->c.dPru->y0;
	}
}

bool L_CalceLista::verifCorrLinDatos(double rMin, double *corrPreg)
{
	L_CalceNodo *cal;
	int i;
	double xx=0,xy=0,yy=0,x=0,y=0;
	double SIGMAxx, SIGMAxy, SIGMAyy, corr;

	for (i=0,cal=root; cal!=NULL; i++,cal=cal->sig)
	{
		xx+=cal->c.dRef->x0*cal->c.dRef->x0;
		xy+=cal->c.dRef->x0*cal->c.dRef->y0;
		yy+=cal->c.dRef->y0*cal->c.dRef->y0;
		x+=cal->c.dRef->x0;
		y+=cal->c.dRef->y0;
	}
	xx/=n;
	xy/=n;
	yy/=n;
	x/=n;
	y/=n;

	SIGMAxx=xx-x*x;
	SIGMAxy=xy-x*y;
	SIGMAyy=yy-y*y;

	if (SIGMAxx>0 && SIGMAyy>0)
		corr=sqrt( (SIGMAxy*SIGMAxy)/(SIGMAxx*SIGMAyy) );
	else
		corr=1;
	if (corrPreg!=NULL)
		*corrPreg=corr;
	if (corr>rMin)
		return false;

	xx=0;
	xy=0;
	yy=0;
	x=0;
	y=0;

	for (i=0,cal=root; cal!=NULL; i++,cal=cal->sig)
	{
		xx+=cal->c.dPru->x0*cal->c.dPru->x0;
		xy+=cal->c.dPru->x0*cal->c.dPru->y0;
		yy+=cal->c.dPru->y0*cal->c.dPru->y0;
		x+=cal->c.dPru->x0;
		y+=cal->c.dPru->y0;
	}
	xx/=n;
	xy/=n;
	yy/=n;
	x/=n;
	y/=n;

	SIGMAxx=xx-x*x;
	SIGMAxy=xy-x*y;
	SIGMAyy=yy-y*y;

	if (SIGMAxx>0 && SIGMAyy>0)
		corr=sqrt( (SIGMAxy*SIGMAxy)/(SIGMAxx*SIGMAyy) );
	else
		corr=1;
	if (corrPreg!=NULL)
		*corrPreg=corr;
	if (corr>rMin)
		return false;
	return true;
}

int L_DescriptorNodo_indLand_qcmp(const void *a, const void *b)
{
	L_DescriptorNodo *da=**(L_DescriptorNodo***)a;
	L_DescriptorNodo *db=**(L_DescriptorNodo***)b;
	if (da->c.indLand > db->c.indLand)
		return 1;
	if (da->c.indLand < db->c.indLand)
		return -1;
	return 0;

}


int L_Calce_id_pru_qcmp(const void *a, const void *b)
{
	const L_Descriptor *da=(*(L_CalceNodo**)a)->c.dPru;
	const L_Descriptor *db=(*(L_CalceNodo**)b)->c.dPru;
	if (da->id > db->id)
		return 1;
	if (da->id < db->id)
		return -1;
	return 0;
}


void L_CalceLista::ordena_id_pru()
{
	sort(0, &L_Calce_id_pru_qcmp);
}

int L_CalceLista::reemplazarReferenciasDescriptores(L_DescriptorLista &desPru)
{
	// desPru ya debe estar ordered por id creciente
	L_CalceNodo *cal;
	L_DescriptorNodo *des;
	int nReem = 0;
	ordena_id_pru();

	for (cal=root,des=desPru.root; cal!=NULL && des!=NULL; cal=cal->sig)
	{
		while (des!=NULL && cal->c.dPru->id > des->c.id)
			des=des->sig;
		if (cal->c.dPru->id == des->c.id)
		{
			cal->c.dPru = &des->c;
			nReem++;
		}
	}
	return nReem;
}


void L_CalceLista::pedirMatrizDatos(L_Matrix &x1y1x2y2)
{
	int i;
	L_CalceNodo *cal;
	x1y1x2y2.reallocate(n, 4);
	for (i=0,cal=root; i<n; i++,cal=cal->sig)
	{
		x1y1x2y2(i,0) = cal->c.dRef->x0;
		x1y1x2y2(i,1) = cal->c.dRef->y0;
		x1y1x2y2(i,2) = cal->c.dPru->x0;
		x1y1x2y2(i,3) = cal->c.dPru->y0;
	}
}

void L_CalceLista::pedirMatrizDatosRayos(L_Matrix &x1y1x2y2, const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru)
{
	int i;
	L_CalceNodo *cal;
	L_RayoCamara r;
	x1y1x2y2.reallocate(n, 4);
	for (i=0,cal=root; i<n; i++,cal=cal->sig)
	{
		cRef.calcRayo(r, cal->c.dRef->x0, cal->c.dRef->y0);
		x1y1x2y2(i,0) = r.tanIzq;
		x1y1x2y2(i,1) = r.tanArr;
		cPru.calcRayo(r, cal->c.dPru->x0, cal->c.dPru->y0);
		x1y1x2y2(i,2) = r.tanIzq;
		x1y1x2y2(i,3) = r.tanArr;
	}
}

void L_CalceLista::pedirMatrizPuntajes(std::vector<double> &puntajes)
{
	L_CalceNodo *cal;
	int i;
	puntajes.resize(size());
	for (i=0, cal=root; i<n; i++, cal=cal->sig)
		puntajes[i] = 1 / cal->c.vectDivDist; // Los descriptores con ratio menor son los mejores
}

void L_CalceLista::eliminaDescriptoresSinLandmark()
{
	L_CalceNodo *tmp=NULL;
	n=0;
	for (pult=&root; *pult!=NULL;)
	{
		if ((*pult)->c.dRef->indLand != -1) // Esta asociado a un landmark
		{
			pult=&(*pult)->sig;
			n++;
		}
		else
		{
			tmp=*pult;
			*pult=(*pult)->sig;
			delete tmp;
		}
	}
}

void L_CalceLista::eliminaLineasHarrisPru(double cornernessMin, double sd, double si, double alfa)
{
	int r = L_GaussianHalfWidth(si);
	int an = 2*r+1;
	int lx, ly, nop_indiv, nop_ant; // Numero de multiplicaciones a hacer con cada opcion
	bool forzarImagenCompleta = false;

	if (root == NULL)
		return;

	lx = root->c.dPru->imorig_lx;
	ly = root->c.dPru->imorig_ly;

	//        convoluciones   gradientes y productos punto a punto
	nop_ant = lx*ly*2*an*4 + lx*ly*5;
	//           convol "sd"   convol "si"  gradientes y productos punto a punto
	nop_indiv = 4*2*an*an*an*n + 2*an*an*n + 4*an*an*5*n;

	if (nop_ant < nop_indiv || forzarImagenCompleta==true)
		eliminaLineasHarrisPru_ant(cornernessMin, sd, si, alfa);
	else
		eliminaLineasHarrisPru_indiv(cornernessMin, sd, si, alfa);
}


void L_CalceLista::eliminaLineasHarrisPru_indiv(double cornernessMin, double sd, double si, double alfa)
{
	if (n <= 0)
		return;
	std::vector<double> x0(n);
	std::vector<double> y0(n);
	L_Array<bool> eliminar(n);
	bool imagen = false;
	int i;
	L_CalceNodo *cal;
	const L_ImageGrayDouble *act = NULL;
	L_ImageGrayDouble imHarris;

	for(i=0, cal=root; i<n && cal!=NULL; i++, cal=cal->sig)
	{
		eliminar[i] = false;
		if(cal->c.dPru->imorig == NULL) // Descriptores artificiales
			continue;
		throw_L_ArgException_if(cal->c.dPru->imorig != act && act!=NULL, "L_CalceLista::eliminaLineasHarrisRef_indiv() : descriptores con distintas imagenes origen");		
		act = cal->c.dPru->imorig;
		if (cal->c.dPru->x0 < 0 || cal->c.dPru->y0 < 0 || cal->c.dPru->x0 > act->lx || cal->c.dPru->y0 > act->ly)
			continue;
		x0[i] = cal->c.dPru->x0;
		y0[i] = cal->c.dPru->y0;
	}
	if (act != NULL) // Por si los descriptores son artificiales
	{
		imHarris.harrisEliminaLineas(cornernessMin, *act, &(x0[0]), &(y0[0]), n, eliminar.data());
		erase_preserving_order(eliminar);
	}
}


void L_CalceLista::eliminaLineasHarrisPru_ant(double cornernessMin, double sd, double si, double alfa)
{
	if (n <= 0)
		return;
	L_Array<bool> eliminar(n);
	int i;
	int x0, y0;
	L_CalceNodo *cal;
	const L_ImageGrayDouble *act = NULL;
	L_ImageGrayDouble imHarris;

	for(i=0, cal=root; i<n && cal!=NULL; i++, cal=cal->sig)
	{
		eliminar[i] = false;
		if(cal->c.dPru->imorig == NULL) // Para descriptores artificiales
			continue;
		if (cal->c.dPru->imorig != act)
		{
			act = cal->c.dPru->imorig;
			imHarris.harris(*act, sd, si, alfa);
		}
		x0 = (int)cal->c.dPru->x0;
		y0 = (int)cal->c.dPru->y0;
		if (x0 < 0)
			x0 = 0;
		if (x0 >= act->lx)
			x0 = act->lx-1;
		if (y0 < 0)
			y0 = 0;
		if (y0 >= act->ly)
			y0 = act->ly-1;
		if (imHarris.pix(x0,y0) < cornernessMin)
			eliminar[i] = true;
	}
	erase_preserving_order(eliminar);
}

int L_CalceLista::eliminaEpipolar_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E)
{
	L_RelEpipolarEsencial rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	int res;
	int i;

	pedirMatrizDatosRayos(xyxy, cRef, cPru);
	pedirMatrizPuntajes(puntajes);
	//res = rel.calcRansac(xyxy, seleccionados, errMax, nHip, nMinAcept, nAceptAutom); // Error: se me olvido poner el errMax... wins
	//res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	//res = rel.calcRansacPreemptivo(xyxy, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	res = rel.calcProsacPreemptivo(xyxy, puntajes, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	E = rel.E;
	return res;
}

int L_CalceLista::eliminaEpipolarPuntajes_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E)
{
	L_RelEpipolarEsencial rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	L_Matrix M, t;
	int res;
	int i;
	L_CalceNodo *cal;

	pedirMatrizDatosRayos(xyxy, cRef, cPru);

	for (i=0,cal=root; i<(int)puntajes.size()&&cal!=NULL; i++,cal=cal->sig)
		puntajes[i] = cal->c.vectDivDist;

	res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	E = rel.E;
	return res;
}

int L_CalceLista::eliminaEpipolarPuntajes_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, L_TransfAfinPI2D &tr, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E)
{
	L_RelEpipolarEsencial rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	L_Matrix M, t;
	double ex, ey;
	int res;
	int i;

	pedirMatrizDatosRayos(xyxy, cRef, cPru);
	tr.calcMatricesRayos(cRef, cPru, M, t);
	
	for (i=0; i<xyxy.li; i++)
	{
		ex = M(0,0) * xyxy(i,0) + M(0,1) * xyxy(i,1) + t(0,0);
		ey = M(1,0) * xyxy(i,2) + M(1,1) * xyxy(i,3) + t(1,0);
		puntajes[i] = - (ex*ex + ey*ey);
	}
	res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	E = rel.E;
	return res;
}

void L_CalceLista::eliminaEpipolar(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, const L_EssentialMatrix &E)
{
	L_Array<bool> elim;
	L_Matrix xyxy;
	L_RelEpipolarEsencial rel;
	int i;
	rel.E = E;
	pedirMatrizDatosRayos(xyxy, cRef, cPru);
	for (i=0; i<n; i++)
	{
		if (rel.errorRelGeom(xyxy(i,0), xyxy(i,1), xyxy(i,2), xyxy(i,3)) > errMax)  // Error: todos decian 0
			elim[i] = true;
		else
			elim[i] = false;
	}
	erase_preserving_order(elim);
}



int L_CalceLista::eliminaEpipolar_Ransac8Puntos(double errMax, int nHip, int nMinAcept, int nAceptAutom, L_MatrizFundamental &F)
{
	L_RelEpipolarFundamental rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	int res;
	int i;

	pedirMatrizDatos(xyxy);
	pedirMatrizPuntajes(puntajes);
	res = rel.calcRansac(xyxy, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	//res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	//res = rel.calcRansacPreemptivo(xyxy, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	//res = rel.calcProsacPreemptivo(xyxy, puntajes, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	F = rel.F;
	return res;
}


int L_CalceLista::eliminaEpipolarPuntajes_Ransac8Puntos(double errMax, int nHip, int nMinAcept, int nAceptAutom, L_MatrizFundamental &F)
{
	L_RelEpipolarFundamental rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	L_Matrix M, t;
	int res;
	int i;
	L_CalceNodo *cal;

	pedirMatrizDatos(xyxy);

	for (i=0,cal=root; i<(int)puntajes.size()&&cal!=NULL; i++,cal=cal->sig)
		puntajes[i] = -cal->c.vectDivDist;
	
	res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	F = rel.F;
	return res;
}


int L_CalceLista::eliminaEpipolarPuntajes_Ransac8Puntos(L_TransfAfinPI2D &tr, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_MatrizFundamental &F)
{
	L_RelEpipolarFundamental rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	L_Matrix M, t;
	double ex, ey;
	int res;
	int i;

	pedirMatrizDatos(xyxy);
	
	for (i=0; i<xyxy.li; i++)
	{
		ex = tr.m11 * xyxy(i,0) + tr.m12 * xyxy(i,1) + tr.tx;
		ey = tr.m21 * xyxy(i,2) + tr.m22 * xyxy(i,3) + tr.ty;
		puntajes[i] = - (ex*ex + ey*ey);
	}

	res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	F = rel.F;
	return res;
}

void L_CalceLista::eliminaEpipolar(double errMax, const L_MatrizFundamental &F)
{
	L_Array<bool> elim;
	L_Matrix xyxy;
	L_RelEpipolarFundamental rel;
	int i;
	rel.F = F;
	pedirMatrizDatos(xyxy);
	for (i=0; i<n; i++)
	{
		if (rel.errorRelGeom(xyxy(i,0), xyxy(i,1), xyxy(i,2), xyxy(i,3)) > errMax) // Error: todos decian 0
			elim[i] = true;
		else
			elim[i] = false;
	}
	erase_preserving_order(elim);
}


int L_CalceLista::eliminaPose_Ransac3Puntos(int proyRef0Pru1, const L_CamaraPinhole &cam, double errMax, int nIter, int nMinAcept, int nAceptAutom, const std::vector<L_CoordsCart3D> &pun, const std::vector<double> &puntajes, L_Pose3D_cuat &pose)
{
	L_TransfPoseProy2D estimadorPose;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	int res;
	int i;

	pedirMatrizDatosRayos(xyxy, cam, cam);

	estimadorPose.ray.resize(xyxy.li);

	if (proyRef0Pru1)
	{
		for (i=0; i<xyxy.li; i++)
		{
			estimadorPose.ray[i].tanIzq = xyxy(i,0);
			estimadorPose.ray[i].tanArr = xyxy(i,1);
		}
	}
	else
	{
		for (i=0; i<xyxy.li; i++)
		{
			estimadorPose.ray[i].tanIzq = xyxy(i,2);
			estimadorPose.ray[i].tanArr = xyxy(i,3);
		}
	}
	estimadorPose.pun = pun;

	res = estimadorPose.calcProsacSimple(puntajes, seleccionados, errMax, nIter, nMinAcept, nAceptAutom);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	pose = estimadorPose.pose;
	return res;
}


int L_CalceLista::eliminaSemejanza_Ransac2puntos(double errMax, int nHip, int nMinAcept, int nAceptAutom, double &rot, double &esc, double &tx, double &ty)
{
	L_TransfSimil2D rel;
	L_Matrix xyxy;
	L_Array<bool> seleccionados(n);
	L_Array<bool> eliminados(n);
	std::vector<double> puntajes(n);
	int res;
	int i;

	pedirMatrizDatos(xyxy);
	pedirMatrizPuntajes(puntajes);
	res = rel.calcRansac(xyxy, seleccionados, errMax, nHip, nMinAcept, nAceptAutom); // Error: se me olvido poner el errMax... wins
	//res = rel.calcProsacSimple(xyxy, puntajes, seleccionados, errMax, nHip, nMinAcept, nAceptAutom); // Error: se me olvido poner el errMax... wins
	//res = rel.calcRansacPreemptivo(xyxy, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	//res = rel.calcProsacPreemptivo(xyxy, puntajes, seleccionados, errMax, nHip, 15*nHip, 1, 3);
	if (res < nMinAcept)
		return res;
	for (i=0; i<n; i++)
		eliminados[i] = !seleccionados[i];
	erase_preserving_order(eliminados);
	rot = rel.ang;
	esc = rel.esc;
	tx = rel.tx;
	ty = rel.ty;
	return res;
}

int L_CalceLista::eliminaDelta(double dMax, double dangMax, double descMax)
{
	n=0;
	const L_Descriptor *d1, *d2;
	L_CalceNodo *ctmp;
	double dist2, dang, desc, dMax2;
	dMax2 = dMax*dMax;
	for (pult=&root; *pult!=NULL; )
	{
		d1 = (*pult)->c.dRef;
		d2 = (*pult)->c.dPru;
		dist2 = (d1->x0 - d2->x0)*(d1->x0 - d2->x0) + (d1->y0 - d2->y0)*(d1->y0 - d2->y0);
		dang = d1->ang - d2->ang;
		if (dang < 0)
			dang = 0;
		desc = d1->sigma0 / d2->sigma0;
		if (desc < 1)
			desc = 1/desc;
		if (dist2 < dMax2 && dang < dangMax && desc < descMax)
		{
			pult = &(*pult)->sig;
			n++;
		}
		else
		{
			ctmp = *pult;
			*pult = (*pult)->sig;
			delete ctmp;
		}
	}
	return n;
}

/*
void L_CalceLista::ejemploRansac(int npuntos)
{
	L_CapturadorImagen capt;
	//L_CapturadorImagen capt("C:/Documents and Settings/Administrador/Mis documentos/Mis vídeos/Bmz/video07.bmz");
	//L_CapturadorImagen capt("C:/Documents and Settings/Administrador/Mis documentos/Mis vídeos/Videos avi slam/aa4.avi");

	L_VentanaImagen vent1(0,0,320,240, "Click para comenzar"), vent2(0, 300, 320,240, "Click para salir"), vent3(0, 600, 320, 240, "todos");
	//L_OpenSURF surf;
#ifdef __L_LOWE_H_
	L_LoweGenDescr surf;
#else
	//L_PDoG_SIFT surf;
#endif
	L_KdTreeBBF kdt;
	L_CalceLista calL;
	L_ImageGrayDouble imEntren, imPrueba;
	L_ImageRGBUchar imCal, imPrincipal;
	L_ShapeArray lins;
	L_EssentialMatrix E;
	L_MatrizFundamental F;
	L_CamaraPinhole cam;
	cam.fijaParametrosPhillips();
	double rot, esc, tx, ty;
	int frame = 0;
	int cant;

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara conectada\n");
		return;
	}
	while(vent1.mouseE() == L_MouseLeido)
	{
		capt.capturarImagen();
		vent1.dibuja(capt.im);
		if (capt.esArchivo())
			break;
		vent1.check();
	}
	vent1.mouseE() = L_MouseLeido;
	
	imEntren = capt.im;
	surf.nobj = 0; // Para la referencia
	surf.calcDescrPuntInt(imEntren);
	surf.desL.eliminaLineasHarris(1.0e-30);
	kdt.createFrom(surf.desL);

	while (vent2.mouseE() == L_MouseLeido)
	{
		if ((frame%15 == 0 && capt.esArchivo()) || vent1.mouseE() != L_MouseLeido)
		{
			vent1.mouseE() = L_MouseLeido;
			surf.desL.clear();
			imEntren = capt.im;
			surf.nobj = 0; // Para la referencia
			surf.calcDescrPuntInt(imEntren);
			surf.desL.eliminaLineasHarris(1.0e-30);
			kdt.destroy();
			kdt.createFrom(surf.desL);
		}
		vent1.check();
		frame++;
		capt.capturarImagen();
		calL.clear();
		surf.desL.clear();
		imPrueba = capt.im;
		surf.nobj = 1; // Para la primera imagen de test
		surf.calcDescrPuntInt(imPrueba);
		//surf.desL.eliminaLineasHarris();
		calL.genCalces(kdt, surf.desL, false, -1);
		cant = (int)calL.size();

		imCal = imEntren;
		calL.dibuja(imCal);
		vent3.dibujaRedimensiona(imCal);

		if (calL.size() < 20)
		{
			printf("CalL: %d -> %d\n", cant, 0);
			continue;
		}
		double t1 = L_TIME();
		switch(npuntos)
		{
		case 0:
			break;
		case 2:
			if (calL.eliminaSemejanza_Ransac2puntos(5.0, 30, 10, (int)(0.9*calL.size()), rot, esc, tx, ty) < 10) // En coord de la imagen
				calL.clear();
			break;
		case 6:
			if (calL.eliminaEpipolar_Ransac6Puntos(cam, cam, 0.005, 32*5*2, 15, (int)(0.9*calL.size()), E) < 20) // En coord normalizadas
				calL.clear();
			break;
		case 8:
			if (calL.eliminaEpipolar_Ransac8Puntos(10, 30, 256, (int)(0.9*calL.size()), F) < 10) // En coord de la imagen
				calL.clear();          // pixels
			break;
		}
		double t2 = L_TIME();
		imPrincipal = capt.im;
		imPrincipal.genDrawing(lins);
		vent1.dibujaRedimensiona(imPrincipal);
		lins.resize(0);
		//printf("CalL: %d -> %d,   %.2f fps\n", cant, calL.size(), 1/(t2-t1));
		imCal = imEntren;
		calL.dibuja(imCal);
		vent2.dibujaRedimensiona(imCal);
		surf.desL.clear();
	}
}

void L_CalceLista::ajustarRansac(int npuntos)
{
	L_CapturadorImagen capt;
	//L_CapturadorImagen capt("C:/Documents and Settings/Administrador/Mis documentos/Mis vídeos/Bmz/video07.bmz");
	//L_CapturadorImagen capt("C:/Documents and Settings/Administrador/Mis documentos/Mis vídeos/Videos avi slam/aa4.avi");

	L_VentanaImagen vent1(0,0,320,200, "Click para comenzar/salir"), ventpf(500, 500, 320,200, "Click para cambiar prefiltrado");
	L_VentanaImagen ventf(500,0,320,200, "Click para cambiar umbral filtrado");

	//L_OpenSURF surf; surf.thres = 0.0001;
	//L_LoweGenDescr surf;
	L_PDoG_SIFT surf;
	L_KdTreeBBF kdt;
	L_CalceLista calL, calLref;
	L_DescriptorLista desLRef, desLpru;
	L_ImageGrayDouble imEntren, imPrueba;
	L_ImageRGBUchar imCal;
	L_ShapeArray lins;
	L_EssentialMatrix E;
	L_MatrizFundamental F;
	L_CamaraPinhole cam;
	cam.fijaParametrosPhillips();
	double rot, esc, tx, ty, errMax = 1;
	int frame = 0;
	int cant;
	bool filtrarCalcesInicio = true;

	if (capt.capturarImagen() == false)
	{
		printf("No hay camara conectada\n");
		return;
	}
	while(vent1.mouseE() == L_MouseLeido)
	{
		capt.capturarImagen();
		vent1.dibuja(capt.im);
		if (capt.esArchivo())
			break;
		vent1.check();
	}
	vent1.mouseE() = L_MouseLeido;
	
	imEntren = capt.im;
	imPrueba = capt.im;
	surf.nobj = 0; // Para la referencia
	surf.calcDescrPuntInt(imEntren);
	surf.desL.eliminaLineasHarris(1.0e-30);
	kdt.createFrom(surf.desL);

	switch(npuntos)
	{
	case 0:
		break;
	case 2:
		errMax = 10;
		break;
	case 6:
		errMax = 0.005;
		break;
	case 8:
		break;
	}

	while (true)
	{
		if ((frame%15 == 0 && capt.esArchivo()) || vent1.mouseE() != L_MouseLeido)
		{
			imEntren = imPrueba;
			vent1.mouseE() = L_MouseLeido;
			surf.desL.clear();
			imEntren = capt.im;
			surf.nobj = 0; // Para la referencia
			surf.calcDescrPuntInt(imEntren);
			surf.desL.eliminaLineasHarris(1.0e-30);
			kdt.destroy();
			kdt.createFrom(surf.desL);
		}

		if (vent1.mouseE() != L_MouseLeido)
			break;

		if (ventf.mouseE() != L_MouseLeido)
		{
			int ex = ventf.mouseX();
			if (ex < 1)
				ex = 1;
			if (ventf.mouseE() == L_MousePresiona)
				errMax *= ex/30.0;
			ventf.mouseE() = L_MouseLeido;
		}

		if (ventpf.mouseE() != L_MouseLeido)
		{
			if (ventpf.mouseE() == L_MousePresiona)
				filtrarCalcesInicio = !filtrarCalcesInicio;
			ventpf.mouseE() = L_MouseLeido;
		}

		vent1.check();
		frame++;
		if (capt.capturarImagen() == false)
			break;
		calL.clear();
		surf.desL.clear();
		vent1.dibujaRedimensiona(capt.im);
		imPrueba = capt.im;
		surf.nobj = 1; // Para la primera imagen de test
		surf.calcDescrPuntInt(imPrueba);
		//surf.desL.eliminaLineasHarris();
		calL.genCalces(kdt, surf.desL, false, -1);
		cant = (int)calL.size();
		if (filtrarCalcesInicio)
			if (calL.eliminaSemejanza_Ransac2puntos(15.0, 30, 10, 40, rot, esc, tx, ty) < 10)
				calL.clear();
		calLref.clear();
		calL.copyListOn(calLref);
		if (calL.size() > 20)
		{
			switch(npuntos)
			{
			case 0:
				break;
			case 2:
				if (calL.eliminaSemejanza_Ransac2puntos(errMax, 30, 10, 40, rot, esc, tx, ty) < 10)
					calL.clear();
				break;
			case 6:
				if (calL.eliminaEpipolar_Ransac6Puntos(cam, cam, errMax, 80, 10, 40, E) < 10)
					calL.clear();
				break;
			case 8:
				if (calL.eliminaEpipolar_Ransac8Puntos(errMax, 30, 15, 40, F) < 15)
					calL.clear();
				break;
			}
		}
		else
			calL.clear();
		printf("CalL: %d -> %d   errLim=%g\n", cant, calL.size(), errMax);
		//
		calLref.dibuja(imCal, imEntren, imPrueba);
		ventpf.dibujaRedimensiona(imCal); // Calces posiblemente prefiltrados
		//
		calL.dibuja(imCal, imEntren, imPrueba);
		lins.resize(0);
		lins.drawLine(30, 0, 30, imCal.ly, 255, 255, 255);
		imCal.genDrawing(lins);
		ventf.dibujaRedimensiona(imCal); // Calces filtrados
		//
		surf.desL.clear();
	}
}
*/


void L_Calcec2Lista::parear(const L_CalceLista &c1, const L_CalceLista &c2)
{
	L_CalceLista v1, v2;
	L_CalceNodo *p1, *p2;
	L_Calcec2 c;
	int cmp;
	c1.copyListOn(v1);
	c2.copyListOn(v2);
	v1.sort(0, &L_Calce::qcmpRef);
	v2.sort(0, &L_Calce::qcmpRef);
	for (p1 = v1.root, p2 = v2.root  ;  p1!=NULL && p2!=NULL; )
	{
		cmp = p1->c.dRef->cmpPosSigmaAng(*p2->c.dRef);
		if (cmp > 0) // Hay que hace avanzar p2, insertar el que avanza
		{
			c.dRef = p2->c.dRef;
			c.dPru = NULL;
			c.dPruc2 = p2->c.dPru;
			push_back(c);
			p2 = p2->sig;
		}
		else if (cmp < 0) // Hay que hacer avanzar p1, insertar el que avanza
		{
			c.dRef = p1->c.dRef;
			c.dPru = p1->c.dPru;
			c.dPruc2 = NULL;
			push_back(c);
			p1 = p1->sig;
		}
		else // Coinciden, insertar ambos pareados y avanzar las 2 listas
		{
			c.dRef = p1->c.dRef;
			c.dPru = p1->c.dPru;
			c.dPruc2 = p2->c.dPru;
			push_back(c);
			p1 = p1->sig;
			p2 = p2->sig;
		}
	}
	// Una lista se acabo, se acabaron los pareados
	// Ahora vienen los de ambas listas separados
	for ( ; p1!=NULL; p1=p1->sig)
	{
		c.dRef = p1->c.dRef;
		c.dPru = p1->c.dPru;
		c.dPruc2 = NULL;
		push_back(c);
	}
	for ( ; p2!=NULL; p2=p2->sig)
	{
		c.dRef = p2->c.dRef;
		c.dPru = NULL;
		c.dPruc2 = p2->c.dPru;
		push_back(c);
	}
	// Se supone que todos los nodos fueron insertados
}

double &L_Transf2D_caract::el(int i)
{
	switch(i)
	{
	case 0: return numVotosTransfFus;
	case 1: return numVotosTransf;
	case 2: return probCelda;
	case 3: return probTransf;
	case 4: return numVotosCelda;
	case 5: return distRomb;
	case 6: return corrLin;
	case 7: return corrPix;
	case 8: return ransac;
	case 9: return calces_1_a_n;
	case 10: return RSLfallados;
	case 11: return porcCerrosCorrectos;
	case 12: return punt;

	default: return numVotosTransfFus;
	}
}

L_String L_Transf2D_caract::imprimeCaract()
{

	return L_String(
		"numVotosTransfFus "
		"numVotosTransf "
		"probCelda "
		"probTransf "
		"numVotosCelda "
		"distRomb "
		"corrLin "
		"corrPix "
		"ransac "
		"calces_1_a_n "
		"RSLfallados "
		"porcCerrosCorrectos "
		"punt"
		);
}

void L_Calcec2Lista::agregar(const L_CalceLista &c1)
{
	L_CalceNodo *p1;
	L_Calcec2 c;
	for (p1 = c1.root ; p1!=NULL; p1=p1->sig)
	{
		c.dRef = p1->c.dRef;
		c.dPru = p1->c.dPru;
		c.dPruc2 = NULL;
		push_back(c);
	}
	// Se supone que todos los nodos fueron insertados
}

void L_TransfPuntInt2D::cuenta_calces_1_a_n()
{
	L_Array<L_CalceNodo *> arr((int)calL.size());
	L_CalceNodo *cal;
	const L_Descriptor *des=NULL;
	int i, calcesConsec=0;
	car.calces_1_a_n=0;
	for (i=0,cal=calL.root ; i<(int)calL.size() ; i++,cal=cal->sig)
		arr[i]=cal;
	qsort(&arr[0], calL.size(), sizeof(arr[0]), L_CalceLista::cmpCalceXY_v2);
	for (i=0; i<(int)calL.size(); i++)
	{
		if (arr[i]->c.dRef!=des)
		{
			des=arr[i]->c.dRef;
			if (calcesConsec>car.calces_1_a_n)
				car.calces_1_a_n=calcesConsec;
			calcesConsec=0;
		}
		calcesConsec++;
	}
	if (calcesConsec>car.calces_1_a_n)
		car.calces_1_a_n=calcesConsec;
}

bool L_TransfPuntInt2D::aplicaRSL(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng)
{
	car.RSLfallados=calL.aplicaRestrSemiLocales(corrigeLista, nVecinos, minCorresp, verifEscalaAngs, varEscala, varAng);
	if (car.RSLfallados>0 && corrigeLista)
		return calcTransfMinCuad();
	return true;
}

bool L_TransfPuntInt2D::aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax)
{
	car.porcCerrosCorrectos=calL.aplicaCuentaCerros(corrigeLista, minDistCerros, maxDistCerros, minCerrosUtil, parecidoCerrosMin, nLineasMax);
	//printf("Cerros correctos: %f\n", porcCerrosCorrectos);
	//if (porcCerrosCorrectos>0 && corrigeLista)
	//	return calcTransfMinCuad();
	return true;
}

int L_TransfPuntInt2D::qcmp(const void *tr1, const void *tr2)
{
	L_TransfPuntInt2D *t1, *t2;
	t1=*(L_TransfPuntInt2D **)tr1;
	t2=*(L_TransfPuntInt2D **)tr2;
	if (t1->car.punt<t2->car.punt)
		return 1;
	if (t1->car.punt>t2->car.punt)
		return -1;
	if (t1->calL.size()<t2->calL.size())
		return 1;
	if (t1->calL.size()>t2->calL.size())
		return -1;
	return 0;
}

int L_TransfPuntInt2DListaPtr::calcNumTransf()
{
	int nMax;
	nMax=0;
	L_TransfPuntInt2DNodoPtr *ptr;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		if (ptr->c==NULL)
			continue;
		nMax++;
	}
	return nMax;
}

int L_TransfPuntInt2DListaPtr::calcNumTransf(L_CountingTree &numTransf)
{
	int n, nMax;
	nMax=0;
	L_TransfPuntInt2DNodoPtr *ptr;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		if (ptr->c==NULL)
			continue;
		else if(ptr->c->calL.root==NULL || ptr->c->calL.root->c.dRef==NULL)
			numTransf.addVoteToFrequencies(-1);
		else
		{
			n=ptr->c->calL.root->c.dRef->nobj;
			numTransf.addVoteToFrequencies(n);
			if (nMax<n)
				nMax=n;
		}
	}
	return nMax;
}

void L_TransfPuntInt2DListaPtr::sort()
{
	L_Array<L_TransfPuntInt2D *> arr(n);
	L_TransfPuntInt2DNodoPtr *ptr;
	int i;
	if (n<2)  // Faltaba esto...
		return;
	for (i=0,ptr=root; i<n; i++,ptr=ptr->sig)
		arr[i]=ptr->c;
	arr.sort(L_TransfPuntInt2D::qcmp);
	for (i=0,ptr=root; i<n; i++,ptr=ptr->sig)
		ptr->c=arr[i];
}

void L_TransfPuntInt2DListaPtr::escribeResumen(FILE *fp)
{
	L_TransfPuntInt2DNodoPtr *ptr;
	L_PuntInt p;
	L_PuntInt p1, p2;
	int lx, ly;
	char tipo[20];

	fprintf(fp,"%d transf\n",n);
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		lx=ptr->c->calL.root->c.dRef->imorig_lx;
		ly=ptr->c->calL.root->c.dRef->imorig_ly;
		p.x0=0;
		p.y0=0;
		ptr->c->proyeccionDe(p, p1);
		p.x0=lx;
		p.y0=ly;
		ptr->c->proyeccionDe(p, p2);
		switch (ptr->c->tipo)
		{
		case L_Tr_Afin:
			strcpy(tipo,"L_Tr_Afin");
			break;
		case L_Tr_Proyectiva:
			strcpy(tipo,"L_Tr_Proyectiva");
			break;
		default:
			strcpy(tipo,"???");
			break;
		}
		fprintf(fp, "vot:%ld punt=%.4f obj=%d (%.1f %.1f)--(%.1f %.1f) R:%d tipo:%s\n", ptr->c->calL.size(), ptr->c->car.punt,
			ptr->c->calL.root->c.dRef->nobj,
			p1.x0, p1.y0, p2.x0, p2.y0,
			(int)ptr->c->car.ransac, tipo);
	}
}

void L_TransfPuntInt2DListaPtr::escribeResumen(const char *nomarch)
{
	FILE *fp;
	fp=fopen(nomarch, "w");
	if (fp==NULL)
		return;
	escribeResumen(fp);
	fclose(fp);
	return;
}

void L_TransfPuntInt2DListaPtr::cuenta_calces_1_a_n()
{
	L_TransfPuntInt2DNodoPtr *ptr;
	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
		ptr->c->cuenta_calces_1_a_n();
}

void L_TransfPuntInt2DListaPtr::aplicaRSL(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng)
{
	L_TransfPuntInt2DNodoPtr *nodo;
	for (pult=&root; *pult!=NULL; )
	{
		if ( (*pult)->c->aplicaRSL(corrigeLista, nVecinos, minCorresp, verifEscalaAngs, varEscala, varAng) == false ) // La transformacion era muy mala
		{
			nodo=*pult;
			*pult=(*pult)->sig;
			delete nodo;
			n--;
		}
		else
			pult=&(*pult)->sig;
	}
}

void L_TransfPuntInt2DListaPtr::aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax, double porcCerrosCorrectosMin)
{
	L_TransfPuntInt2DNodoPtr *nodo;
	for (pult=&root; *pult!=NULL; )
	{
		(*pult)->c->aplicaCuentaCerros(corrigeLista, minDistCerros, maxDistCerros, minCerrosUtil, parecidoCerrosMin, nLineasMax);
		if ( corrigeLista && (*pult)->c->car.porcCerrosCorrectos <  porcCerrosCorrectosMin ) // La transformacion era muy mala
		{
			nodo=*pult; // Eliminar la transformacion
			*pult=(*pult)->sig;
			delete nodo;
			n--;
		}
		else
			pult=&(*pult)->sig;
	}
}

// Aca en lins el numero de imagen es: 0 = test,  1 = referencia 0,  2 = referencia 1,  etc.
void L_TransfPuntInt2DListaPtr::dibujaTransformaciones(L_ShapeArray &lins, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica, int transfEspecifica)
{
	L_CountingTree numTransf;
	L_TransfPuntInt2DNodoPtr *ptr;
	L_PuntInt pRef1, pPru1, pRef2, pPru2;
	int nImagen, nTransf, nImagenTot;
	int lx, ly;
	int i, j;
	L_uchar R,G,B;

	for (ptr=root; ptr!=NULL; ptr=ptr->sig)
	{
		if (ptr->c->calL.root==NULL || ptr->c->calL.root->c.dRef==NULL) // No deberia pasar, pero ...
		{
			printf("Situacion rara en L_TransfPuntInt2DListaPtr::dibujaTransformaciones:\n");
			printf("Transformacion sin lista de calces (asi no se puede dibujar)\n");
			continue;
		}
		// Se tiene la transformacion ptr->cont. Hay que mostrar las mejores transformaciones que cumplan con el criterio

		// Si maxImRef es negativo, la cantidad de imagenes de referencia a mostrar es libre. Si no, es acotada
		// Si maxTransfPorImagen es negativo, la cantidad de transformaciones por imagen de referencia a mostrar es libre. Si no, es acotada
		// Si maxTransfTotales es negativo, la cantidad de transformaciones totales a mostrar es libre. Si no, es acotada
		// Si imRefEspecifica es mayor o igual a zero, hay que mostrar esa imagen en particular
		// Si transfEspecifica es mayor o igual a zero, hay que mostrar ese numero de transformacion en particular

		if (maxTransfTotales>=0 && numTransf.nVotTot>=maxTransfTotales) // Se tienen todas las transformaciones pedidas
			break;

		nImagen=ptr->c->calL.root->c.dRef->nobj; // Numero de imagen actual
		nTransf=numTransf.votosDe(nImagen);       // Numero de transformaciones agregadas para esta imagen
		nImagenTot=numTransf.nCellsTot;  // Numero de imagenes ya agregadas a la lista

		if (imRefEspecifica>=0 && nImagen!=imRefEspecifica) // Se pide una imagen específica
			continue;

		if (transfEspecifica>=0 && nTransf<transfEspecifica) // Falta para la transformacion pedida
			continue;
		if (transfEspecifica>=0 && nTransf>transfEspecifica) // Ya se paso por la transformacion pedida
			break;

		if (maxImRef>=0 && nImagenTot>=maxImRef && nTransf==0) // No se pueden agregar nuevas imagenes a la lista
			continue;
		if (maxTransfPorImagen>=0 && nTransf>=maxTransfPorImagen) // Para esta imagen ya hay muchas transformaciones
			continue;

		// Paso todos los filtros, hay que agregarla
		if (calces1_o_caja0==true)
			ptr->c->calL.genLineas(lins, nImagen+1); // Es una lata sumarle 1 pero es porque imPru = 0, imRef1 = 1, etc...
		else
		{
			lx=ptr->c->calL.root->c.dRef->imorig_lx;
			ly=ptr->c->calL.root->c.dRef->imorig_ly;
			lins.genRandomColor(R, G, B);
			pRef1.x0=0;  pRef1.y0=0;  ptr->c->proyeccionDe(pRef1, pPru1);
			pRef2.x0=lx; pRef2.y0=0;  ptr->c->proyeccionDe(pRef2, pPru2);
			// Dibujar 1 linea gruesa = 9 lineas finas
			for (i=-1; i<=1; i++)
				for (j=-1; j<=1; j++)
					lins._drawShape(0,0,(int)(pPru1.x0+0.5)+i, (int)(pPru1.y0+0.5)+j, (int)(pPru2.x0+0.5)+i, (int)(pPru2.y0+0.5)+j, R, G, B);
			pRef1.x0=lx; pRef1.y0=ly; ptr->c->proyeccionDe(pRef1, pPru1);
			lins._drawShape(0,0,(int)(pPru1.x0+0.5), (int)(pPru1.y0+0.5), (int)(pPru2.x0+0.5), (int)(pPru2.y0+0.5), R, G, B);
			pRef2.x0=0;  pRef2.y0=ly; ptr->c->proyeccionDe(pRef2, pPru2);
			lins._drawShape(0,0,(int)(pPru1.x0+0.5), (int)(pPru1.y0+0.5), (int)(pPru2.x0+0.5), (int)(pPru2.y0+0.5), R, G, B);
			pRef1.x0=0;  pRef1.y0=0;  ptr->c->proyeccionDe(pRef1, pPru1);
			lins._drawShape(0,0,(int)(pPru1.x0+0.5), (int)(pPru1.y0+0.5), (int)(pPru2.x0+0.5), (int)(pPru2.y0+0.5), R, G, B);
			pRef2.x0=lx/2, pRef2.y0=ly/2; ptr->c->proyeccionDe(pRef2, pPru2);
			lins._drawShape(nImagen+1,0,(int)(pRef2.x0+0.5),(int)(pRef2.y0+0.5),(int)(pPru2.x0+0.5),(int)(pPru2.y0+0.5),R,G,B);
		}
		numTransf.addVoteToFrequencies(nImagen);
	}
}

L_TransfPuntInt2D *L_TransfPuntInt2DListaPtr::buscaTransf(int nobj, int indice)
{
	L_TransfPuntInt2DNodoPtr *nodo;
	int i = 0;
	for (nodo=root; nodo!=NULL; nodo=nodo->sig)
	{
		if (nodo->c==NULL)
		{
			printf("Error en L_TransfPuntInt2DListaPtr::buscaTransf(): lista mal formada\n");
			return NULL;
		}
		if (nodo->c->pideNumObj() == nobj)
		{
			if (i == indice)
				return nodo->c;
			else
				i++;
		}
	}
	return NULL;
}

void L_TransfPuntInt2DListaPtr::destroyList(L_List<L_Calce> *cal_pool)
{
    L_NodePtr<L_TransfPuntInt2D> *p;
	while(root!=NULL)
	{
		p=root->sig;
		if (cal_pool != NULL)
			root->c->calL.moveListTo(*cal_pool);
		delete root->c;
		delete root;
		root=p;
	}
	pult=&root;
	n=0;
}

bool L_TransfAfinPI2D::calcTransfExacto(const L_Calce calArr[])
{
	L_Matrix A,x,b;
	int i;
	A.reallocate(6,6);
	b.reallocate(6,1);
	for (i=0; i<3; i++)
	{
		A(2*i+0,0)=calArr[i].dRef->x0;
		A(2*i+0,1)=calArr[i].dRef->y0;
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=calArr[i].dRef->x0;
		A(2*i+1,4)=calArr[i].dRef->y0;
		A(2*i+1,5)=1;

		b(2*i+0,0)=calArr[i].dPru->x0;
		b(2*i+1,0)=calArr[i].dPru->y0;
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

bool L_TransfAfinPI2D::calcTransfMinCuad()
{
	L_PUSH_EXECUTING_FN("L_TransfAfinPI2D::calcTransfMinCuad");
	L_Matrix A,AT,ATA,ATb,x,b;
	L_CalceNodo *nod;
	int i;

	A.reallocate(2*(int)calL.size(),6);
	b.reallocate(2*(int)calL.size(),1);
	nod=calL.root;
	for (i=0; i<(int)calL.size(); i++)
	{
		throw_L_ArgException_if(nod==NULL,"L_TransfAfinPI2D::calcTransfMinCuad");
		A(2*i+0,0)=nod->c.dRef->x0;
		A(2*i+0,1)=nod->c.dRef->y0;
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=nod->c.dRef->x0;
		A(2*i+1,4)=nod->c.dRef->y0;
		A(2*i+1,5)=1;

		b(2*i+0,0)=nod->c.dPru->x0;
		b(2*i+1,0)=nod->c.dPru->y0;
		nod=nod->sig;
	}
	AT.transpOf(A);
	ATA.OP_mult(AT,A);
	ATb.OP_mult(AT,b);

	ATb.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(ATA) == false)
	{
		L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcTransfMinCuad");
		return false;
	}

	m11=x(0,0);
	m12=x(1,0);
	tx=x(2,0);
	m21=x(3,0);
	m22=x(4,0);
	ty=x(5,0);
	L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcTransfMinCuad");
	return true;
}

bool L_TransfAfinPI2D::calcParametrosInversosDe(const L_TransfAfinPI2D &other)
{
	L_Matrix A;
	A.reallocate(3,3);
	A(0,0)=other.m11;
	A(0,1)=other.m12;
	A(0,2)=other.tx;
	A(1,0)=other.m21;
	A(1,1)=other.m22;
	A(1,2)=other.ty;
	A(2,0)=0;
	A(2,1)=0;
	A(2,2)=1;
	if (A.invertMe()==false)
		return false;
	m11=A(0,0);
	m12=A(0,1);
	tx=A(0,2);
	m21=A(1,0);
	m22=A(1,1);
	ty=A(1,2);
	return true;
}

L_TransfAfinPI2D *L_TransfAfinPI2D::clona()
{
	L_TransfAfinPI2D *a=new L_TransfAfinPI2D;
	if (a==NULL) return NULL;
	a->m11=m11;
	a->m12=m12;
	a->tx=tx;
	a->m21=m21;
	a->m22=m22;
	a->ty=ty;
	calL.copyListOn(a->calL);
	a->L_Transf2D_POD::operator=(*this);

	return a;
}

L_TransfAfinPI2D *L_TransfAfinPI2D::clonaRegalaCalces()
{
	L_TransfAfinPI2D *a=new L_TransfAfinPI2D;
	if (a==NULL) return NULL;
	a->m11=m11;
	a->m12=m12;
	a->tx=tx;
	a->m21=m21;
	a->m22=m22;
	a->ty=ty;
	calL.moveListTo(a->calL);
	a->L_Transf2D_POD::operator=(*this);

	return a;
}

void L_TransfAfinPI2D::swap(L_TransfAfinPI2D &other)
{
	L_Transf2D_POD t2d_t;
	L_TransfAfin2D ta_t;

	t2d_t.L_Transf2D_POD::operator=(other);
	other.L_Transf2D_POD::operator=(*this);
	this->L_Transf2D_POD::operator=(t2d_t);

	ta_t.L_TransfAfin2D::operator=(other);
	other.L_TransfAfin2D::operator=(*this);
	this->L_TransfAfin2D::operator=(ta_t);

	other.calL.swap(calL);
}

bool L_TransfAfinPI2D::calcTransfRansac(double errX, double errY, int nMinAcept, int nIntentos)
{
	L_PUSH_EXECUTING_FN("L_TransfAfinPI2D::calcTransfRansac");
	L_Array<bool> compat((int)calL.size());
	L_Array<L_CalceNodo*> arr((int)calL.size());
	L_TransfAfinPI2D invTr;
	int i, j, k;
	int azar, nacept;
	L_CalceNodo *cal;
	double err;

#ifdef L_BASIC_DEBUG
	if (calL.size()<=0)
		L_hard_shutdown("L_TransfAfinPI2D::calcTransfRansac : non-positive size allocation");
#endif

	for (i=0, cal=calL.root; i<(int)calL.size(); i++,cal=cal->sig)
		arr[i]=cal;

	for (i=0; i<nIntentos; i++)
	{
		// Resetear los calces compatibles
		for (j=0; j<(int)calL.size(); j++)
			compat[j]=false;
		nacept=0;
		// Seleccionar calces iniciales
		k=0;
		for (j=0; j<3; j++)
		{
			azar=(int)( rand()%calL.size() );
			if (compat[azar]==true) // mala suerte
			{
				j--;
				k++;
				continue;
			}
			if (k>=10)
			{
				k=(i*200+i*i*10) % ((int)calL.size()-3);
				compat[k]=true;
				compat[k+1]=true;
				compat[k+2]=true;
				break;
			}
			compat[azar]=true;
			nacept++;
			//invTr.calL.insertaNodoInverso(*arr[i]); // No puedo creer que eso haya estado aca !!!!!
			invTr.calL.insertaNodoInverso(arr[azar]->c);
		}
		if (invTr.calcTransfMinCuad()==false)
			continue;
		for (j=0; j<(int)calL.size(); j++)
		{
			if (compat[j]==true)
				continue;
			// errores de retroproyeccion
			err=arr[j]->c.dRef->x0 - (invTr.m11*arr[j]->c.dPru->x0 + invTr.m12*arr[j]->c.dPru->y0 + invTr.tx);
			if (err<0)
			{
				if (-err > errX)
					continue;
			}
			else
			{
				if (err > errX)
					continue;
			}
			err=arr[j]->c.dRef->y0 - (invTr.m21*arr[j]->c.dPru->x0 + invTr.m22*arr[j]->c.dPru->y0+ invTr.ty);
			if (err<0)
			{
				if (-err > errY)
					continue;
			}
			else
			{
				if (err > errY)
					continue;
			}
			compat[j]=true;
			nacept++;
		}
		if (nacept>=nMinAcept)
		{
			calL.root=NULL;
			calL.pult=&calL.root;
			for (j=0; j<(int)calL.size(); j++)
			{
				if (compat[j]==true)
				{
					*calL.pult=arr[j];
					calL.pult=&(*calL.pult)->sig;
				}
				else
					delete(arr[j]);
			}
			*calL.pult=NULL;
			calL.n=nacept;
			calcTransfMinCuad();
			L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcTransfRansac");
			return true;
		}
		invTr.calL.clear();
	}

	L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcTransfRansac");
	return false;
}

bool L_TransfAfinPI2D::deformaImagen(L_ImageGrayDouble &dest, const L_ImageGrayDouble &orig, L_Transf2DDeformTipo tipo) const
{
	L_PuntInt p1, p2;
	L_TransfAfinPI2D itr;
	double x=0, y=0;
	int i;
	double punBorde[4][2]={ {0,0} , {1,0}, {1,1}, {0,1} };
	double nPunBorde=4;

	switch(tipo)
	{
	case PIXELADO:
	case BILIN_1NIV_ISOTR:
	case BILIN_XNIV_ISOTR:
		break;
	case BILIN_1NIV_ANISOTR:
	case BILIN_XNIV_ANISOTR:
	default:
		return false;
	}
	for (i=0; i<nPunBorde; i++)
	{
		p1.x0 = punBorde[i][0]*orig.lx;
		p1.y0 = punBorde[i][1]*orig.ly;
		proyeccionDe(p1, p2);
		if (p2.x0>x)
			x=p2.x0;
		if (p2.y0>y)
			y=p2.y0;
	}
	dest.reallocate((int)x, (int)y);
	if (!itr.calcParametrosInversosDe(*this))
		return false;
	switch(tipo)
	{
	case PIXELADO:
		dest.AplicaTransfAfinPixelado(orig, this->m11, this->m12, this->m21, this->m22, this->tx, this->ty, false);
		break;

	default:
		return false; // Hay que implementarlo
	}
	return true;
}

bool L_TransfAfinPI2D::calcVarParametros(L_Matrix &var)
{
	L_PUSH_EXECUTING_FN("L_TransfAfinPI2D::calcSigmaParametros");
	L_Matrix A,AT,ATA,ATb,x,b;
	L_CalceNodo *nod;
	int i;
	double xPr, yPr, e=0;
	A.reallocate(2*(int)calL.size(),6);

	b.reallocate(2*(int)calL.size(),1);
	nod=calL.root;
	for (i=0; i<(int)calL.size(); i++)
	{
		throw_L_ArgException_if(nod==NULL, "L_TransfAfinPI2D::calcVarParametros");
		A(2*i+0,0)=nod->c.dRef->x0;
		A(2*i+0,1)=nod->c.dRef->y0;
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=nod->c.dRef->x0;
		A(2*i+1,4)=nod->c.dRef->y0;
		A(2*i+1,5)=1;

		b(2*i+0,0)=nod->c.dPru->x0;
		b(2*i+1,0)=nod->c.dPru->y0;
		nod=nod->sig;
	}
	AT.transpOf(A);
	ATA.OP_mult(AT,A);
	if (ATA.invertMe()==false)
	{
		L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcTransfMinCuad");
		return false;
	}
	ATb.OP_mult(AT,b);
	x.OP_mult(ATA,ATb);
	m11=x(0,0);
	m12=x(1,0);
	tx=x(2,0);
	m21=x(3,0);
	m22=x(4,0);
	ty=x(5,0);
	// Calculo de sigmas
	nod=calL.root;
	for (i=0; i<(int)calL.size(); i++)
	{
		throw_L_ArgException_if(nod==NULL, "L_TransfAfinPI2D::calcVarParametros");
		xPr=this->m11*nod->c.dRef->x0 + this->m12*nod->c.dRef->y0 + this->tx;
		yPr=this->m21*nod->c.dRef->x0 + this->m22*nod->c.dRef->y0 + this->ty;
		e+=(xPr-nod->c.dPru->x0)*(xPr-nod->c.dPru->x0) + (yPr-nod->c.dPru->y0)*(yPr-nod->c.dPru->y0);
		nod=nod->sig;
	}
	e/=(calL.size()-3); // Se requieren 3 calces para fijar la transf -> n-3 grados de libertad
	ATA*=e;
	ATA.swap(var);
	L_POP_EXECUTING_FN("L_TransfAfinPI2D::calcSigmaParametros");
	return true;
}

void L_TransfAfinPI2D::linealiza(L_TransfPuntInt2D *tr, double x0, double y0)
{
	L_PuntInt o0, o1, o2;
	L_PuntInt t0, t1, t2;
	double d = 1.0e-5;
	o0.x0 = x0;   o0.y0 = y0;
	o1.x0 = x0+d; o1.y0 = y0;
	o2.x0 = x0;   o2.y0 = y0+d;
	tr->proyeccionDe(o0,t0); // t0_x = m11*  x0   + m12*  y0   + tx,   t0_y = m21*  x0   + m22*  y0   + ty
	tr->proyeccionDe(o1,t1); // t1_x = m11*(x0+d) + m12*  y0   + tx,   t1_y = m21*(x0+d) + m22*  y0   + ty
	tr->proyeccionDe(o2,t2); // t2_x = m11*  x0   + m12*(y0+d) + tx,   t2_y = m21*  x0   + m22*(y0+d) + ty

	tx = (d*t0.x0-x0*t1.x0+x0*t0.x0-y0*t2.x0+y0*t0.x0)/d;
	ty = (d*t0.y0-x0*t1.y0+x0*t0.y0-y0*t2.y0+y0*t0.y0)/d;
	m11 = -(-t1.x0+t0.x0)/d;
	m12 = -(-t2.x0+t0.x0)/d;
	m21 = -(-t1.y0+t0.y0)/d;
	m22 = -(-t2.y0+t0.y0)/d;
}

/*
void L_TransfAfinPI2D::descomponeSVD(double &escMax, double &escMin, double &rot, const L_CamaraPinhole *cR, const L_CamaraPinhole *cP)
{
	L_Matrix M(2,2), t(2,1), U(2,2),D(2,2),V(2,2), VT(2,2), R(2,2);
	double ang1, ang2;
	if (cR != NULL && cP!=NULL)
		calcMatricesRayos(*cR, *cP, M, t);
	else
	{
		M[0][0] = m11;
		M[0][1] = m12;
		M[1][0] = m21;
		M[1][1] = m22;
		t[0][0] = tx;
		t[1][0] = ty;
	}
	U[0][0] = M[0][0];
	U[0][1] = M[0][1];
	U[1][0] = M[1][0];
	U[1][1] = M[1][1];
	if (M.svd2x2(ang1, ang2, escMax, escMin)
	{
		rot = ang1 + ang2;
		if (escMin > escMax)
		{
			double temp = escMax;
			escMax = escMin;
			escMin = temp;
		}
	}
	else
	{
		U.svd_ordenado(D,V);
		VT.transpOf(V);
		R.OP_mult(U,VT);
		escMax = D[0][0];
		escMin = D[1][0];
		rot = atan2(R[0][1], R[0][0]);
	}
}
*/


bool L_TransfProyectivaPI2D::calcTransfExacto(const L_Calce *calArr)
{
#if defined(x) || defined(y) || defined(u) || defined(v)
#error A quien se le ocurre dejar definido x, y, u, v .... .. .  .
#endif
	L_Matrix A, x, b;
	int i;
	A.reallocate(4*2,8);
	b.reallocate(4*2,1);
	for (i=0; i<4; i++)
	{
#define x ( calArr[i].dRef->x0 )  // Va de (x,y) a (u,v)
#define y ( calArr[i].dRef->y0 )
#define u ( calArr[i].dPru->x0 )
#define v ( calArr[i].dPru->y0 )
		A(2*i,0) = x;
		A(2*i,1) = y;
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = -x*u;
		A(2*i,7) = -y*u;
		b(2*i,0) = u;
		//
		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = x;
		A(2*i+1,4) = y;
		A(2*i+1,5) = 1;
		A(2*i+1,6) = -x*v;
		A(2*i+1,7) = -y*v;
		b(2*i+1,0) = v;
#undef x
#undef y
#undef u
#undef v
	}

	b.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(A) == false)
		return false;

	for (i=0; i<8; i++)
		el(i)=x(i,0);
	return true;
}

bool L_TransfProyectivaPI2D::calcTransfExacto(double *xR, double *yR, double *xP, double *yP)
{
#if defined(x) || defined(y) || defined(u) || defined(v)
#error A quien se le ocurre dejar definido x, y, u, v .... .. .  .
#endif
	L_Matrix A, x, b;
	int i;
	A.reallocate(4*2,8);
	b.reallocate(4*2,1);
	for (i=0; i<4; i++)
	{
#define x ( xR[i] )  // Va de (x,y) a (u,v)
#define y ( yR[i] )
#define u ( xP[i] )
#define v ( yP[i] )
		A(2*i,0) = x;
		A(2*i,1) = y;
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = -x*u;
		A(2*i,7) = -y*u;
		b(2*i,0) = u;
		//
		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = x;
		A(2*i+1,4) = y;
		A(2*i+1,5) = 1;
		A(2*i+1,6) = -x*v;
		A(2*i+1,7) = -y*v;
		b(2*i+1,0) = v;
#undef x
#undef y
#undef u
#undef v
	}
/*
	if (A.invertMe() == false)
		return false;
	x.OP_mult(A,b);
*/
	b.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(A) == false)
		return false;

	for (i=0; i<8; i++)
		el(i)=x(i,0);
	return true;
}

bool L_TransfProyectivaPI2D::calcTransfMinCuad()
{
#if defined(x) || defined(y) || defined(u) || defined(v)
#error A quien se le ocurre dejar definido x, y, u, v .... .. .  .
#endif
	L_Matrix A, x, b;
	L_Matrix AT, ATA, ATb;
	L_CalceNodo *c;
	int i;
	A.reallocate((int)calL.size()*2,8);
	b.reallocate((int)calL.size()*2,1);
	for (i=0, c=calL.root; i<(int)calL.size() && c!=NULL; i++, c=c->sig)
	{
#define x ( c->c.dRef->x0 )  // Va de (x,y) a (u,v)
#define y ( c->c.dRef->y0 )
#define u ( c->c.dPru->x0 )
#define v ( c->c.dPru->y0 )
		A(2*i,0) = x;
		A(2*i,1) = y;
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = -x*u;
		A(2*i,7) = -y*u;
		b(2*i,0) = u;
		//
		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = x;
		A(2*i+1,4) = y;
		A(2*i+1,5) = 1;
		A(2*i+1,6) = -x*v;
		A(2*i+1,7) = -y*v;
		b(2*i+1,0) = v;
#undef x
#undef y
#undef u
#undef v
	}
	AT.transpOf(A);
	ATb.OP_mult(AT, b);
	ATA.OP_mult(AT, A);

	ATb.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(ATA) == false)
		return false;

	for (i=0; i<8; i++)
		el(i)=x(i,0);
	return true;
}

bool L_TransfProyectivaPI2D::calcTransfMinCuad(double *xR, double *yR, double *xP, double *yP, int n)
{
#if defined(x) || defined(y) || defined(u) || defined(v)
#error A quien se le ocurre dejar definido x, y, u, v .... .. .  .
#endif
	L_Matrix A, x, b;
	L_Matrix AT, ATA, ATb;
	int i;
	A.reallocate(n*2,8);
	b.reallocate(n*2,1);
	for (i=0; i<(int)calL.size(); i++)
	{
#define x ( xR[i] )  // Va de (x,y) a (u,v)
#define y ( yR[i] )
#define u ( xP[i] )
#define v ( yP[i] )
		A(2*i,0) = x;
		A(2*i,1) = y;
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = -x*u;
		A(2*i,7) = -y*u;
		b(2*i,0) = u;
		//
		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = x;
		A(2*i+1,4) = y;
		A(2*i+1,5) = 1;
		A(2*i+1,6) = -x*v;
		A(2*i+1,7) = -y*v;
		b(2*i+1,0) = v;
#undef x
#undef y
#undef u
#undef v
	}
	AT.transpOf(A);
	ATb.OP_mult(AT, b);
	ATA.OP_mult(AT, A);

	ATb.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(ATA) == false)
		return false;

	for (i=0; i<8; i++)
		el(i)=x(i,0);
	return true;
}

bool L_TransfProyectivaPI2D::calcParametrosInversosDe(const L_TransfProyectivaPI2D &other)
{
	L_Matrix M;
	int i;
	M.reallocate(3,3);
	for (i=0; i<8; i++)
		M(i/3,i%3)=other.el(i);
	M(2,2)=1;
	if (M.invertMe() == false)
		return false;
	for (i=0; i<8; i++)
		el(i)=M(i/3,i%3)/M(2,2);
	return true;
}


bool L_TransfProyectivaPI2D::calcTransfRansac_2cal_ant(double errX, double errY, int nMinAcept, int nIntentos)
{
	L_TransfProyectivaPI2D invTr;
	L_Array<bool> compat((int)calL.size());
	L_Array<L_CalceNodo*> arr((int)calL.size());
	L_CalceNodo *cal;
	L_Array<L_CalceNodo> arrTmp((int)calL.size());
	std::vector<L_Descriptor> desTmp((int)calL.size());
	double den, err;
	int azar, nacept;
	int i, j;

#ifdef L_BASIC_DEBUG
	if (calL.size()<=0)
		L_hard_shutdown("L_TransfProyectivaPI2D::calcTransfRansac : non-positive size allocation");
#endif
	if ((int)calL.size() < nMinAcept) // Asi ni modo...
		return false;

	for (i=0, cal=calL.root; i<(int)calL.size(); i++,cal=cal->sig)
		arr[i]=cal;

	// Generar calces artificiales extra
	for (i=0; i<(int)calL.size(); i++)
	{
		desTmp[i*2].x0 = arrTmp[i-(int)calL.size()].c.dPru->x0PuntaFlecha();  // Pru x
		desTmp[i*2].y0 = arrTmp[i-(int)calL.size()].c.dPru->y0PuntaFlecha();  // Pru y
		desTmp[i*2+1].x0 = arrTmp[i-(int)calL.size()].c.dRef->x0PuntaFlecha(); // Ref x
		desTmp[i*2+1].y0 = arrTmp[i-(int)calL.size()].c.dRef->y0PuntaFlecha(); // Ref y
		arrTmp[i].c.dPru=&desTmp[i*2];
		arrTmp[i].c.dRef=&desTmp[i*2+1];
		arr[i+(int)calL.size()]=&arrTmp[i];
	}

	for (i=0; i<nIntentos; i++)
	{
		// Resetear los calces compatibles
		for (j=0; j<(int)calL.size(); j++)
			compat[j]=false;
		// Seleccionar calces iniciales para generar la transformacion (incluyendo los artificiales)
		nacept=0;
		while (nacept < 2)
		{
			azar=(int)( rand()%calL.size() ); // Hay que mejorar esto !!!!
			if (compat[azar]==true) // mala suerte
				continue;
			compat[azar]=true;
			nacept+=1;
			invTr.calL.insertaNodoInverso(arr[azar]->c);
			invTr.calL.insertaNodoInverso(arr[azar+(int)calL.size()]->c); // Esta es la unica linea en que se usan los calces artificiales
		}
		if (invTr.calcTransfMinCuad()==false)
			continue;
		// Los calces artificiales no cuentan para el consenso...
		for (j=0; j<(int)calL.size(); j++)
		{
			if (compat[j]==true)
				continue;
			// errores de retroproyeccion
			den = invTr.m31*arr[j]->c.dPru->x0 + invTr.m32*arr[j]->c.dPru->y0 + 1; //m33 = 1
			err=arr[j]->c.dRef->x0 - (invTr.m11*arr[j]->c.dPru->x0 + invTr.m12*arr[j]->c.dPru->y0 + invTr.tx)/den;
			if (err<0)
			{
				if (-err > errX)
					continue;
			}
			else
			{
				if (err > errX)
					continue;
			}
			err=arr[j]->c.dRef->y0 - (invTr.m21*arr[j]->c.dPru->x0 + invTr.m22*arr[j]->c.dPru->y0 + invTr.ty)/den;
			if (err<0)
			{
				if (-err > errY)
					continue;
			}
			else
			{
				if (err > errY)
					continue;
			}
			compat[j]=true;
			nacept++;
		}
		if (nacept>=nMinAcept)
		{
			calL.root=NULL;
			calL.pult=&calL.root;
			for (j=0; j<(int)calL.size(); j++)
			{
				if (compat[j]==true)
				{
					*calL.pult=arr[j];
					calL.pult=&(*calL.pult)->sig;
				}
				else
					delete(arr[j]);
			}
			*calL.pult=NULL;
			calL.n=nacept;
			calcTransfMinCuad();
			return true;
		}
		invTr.calL.clear();
	}
	return false;
}

bool L_TransfProyectivaPI2D::calcTransfRansac_4cal_ant(double errX, double errY, int nMinAcept, int nIntentos)
{
	L_TransfProyectivaPI2D invTr;
	L_Array<bool> compat((int)calL.size());
	int i, j;
	int azar, nacept;
	L_Array<L_CalceNodo*> arr((int)calL.size());
	L_CalceNodo *cal;
	double den, err;

#ifdef L_BASIC_DEBUG
	if (calL.size()<=0)
		L_hard_shutdown("L_TransfProyectivaPI2D::calcTransfRansac : non-positive size allocation");
#endif
	if ((int)calL.size() < nMinAcept)
		return false;

	for (i=0, cal=calL.root; i<(int)calL.size(); i++,cal=cal->sig)
		arr[i]=cal;

	for (i=0; i<nIntentos; i++)
	{
		// Resetear los calces compatibles
		for (j=0; j<(int)calL.size(); j++)
			compat[j]=false;
		nacept=0;
		// Seleccionar calces iniciales
		nacept=0;
		while (nacept < 4)
		{
			azar=(int)( rand()%calL.size() );
			if (compat[azar]==true) // mala suerte
				continue;
			compat[azar]=true;
			nacept++;
			invTr.calL.insertaNodoInverso(arr[azar]->c);
		}
		if (invTr.calcTransfMinCuad()==false)
			continue;
		for (j=0; j<(int)calL.size(); j++)
		{
			if (compat[j]==true)
				continue;
			// errores de retroproyeccion
			den = invTr.m31*arr[j]->c.dPru->x0 + invTr.m32*arr[j]->c.dPru->y0 + 1; //m33 = 1
			err=arr[j]->c.dRef->x0 - (invTr.m11*arr[j]->c.dPru->x0 + invTr.m12*arr[j]->c.dPru->y0 + invTr.tx)/den;
			if (err<0)
			{
				if (-err > errX)
					continue;
			}
			else
			{
				if (err > errX)
					continue;
			}
			err=arr[j]->c.dRef->y0 - (invTr.m21*arr[j]->c.dPru->x0 + invTr.m22*arr[j]->c.dPru->y0 + invTr.ty)/den;
			if (err<0)
			{
				if (-err > errY)
					continue;
			}
			else
			{
				if (err > errY)
					continue;
			}
			compat[j]=true;
			nacept++;
		}
		if (nacept>=nMinAcept)
		{
			calL.root=NULL;
			calL.pult=&calL.root;
			for (j=0; j<(int)calL.size(); j++)
			{
				if (compat[j]==true)
				{
					*calL.pult=arr[j];
					calL.pult=&(*calL.pult)->sig;
				}
				else
					delete(arr[j]);
			}
			*calL.pult=NULL;
			calL.n=nacept;
			calcTransfMinCuad();
			return true;
		}
		invTr.calL.clear();
	}
	return false;
}

bool L_TransfProyectivaPI2D::calcTransfRansac_ant(double errX, double errY, int nMinAcept, int nIntentos)
{
	bool ret;
	L_PUSH_EXECUTING_FN("L_TransfProyectivaPI2D::calcTransfRansac");
	if (ransac_2cal == true)
		ret=calcTransfRansac_2cal_ant(errX, errY, nMinAcept, nIntentos);
	else
		ret=calcTransfRansac_4cal_ant(errX, errY, nMinAcept, nIntentos);
	L_POP_EXECUTING_FN("L_TransfProyectivaPI2D::calcTransfRansac");
	return ret;
}

bool L_TransfProyectivaPI2D::calcTransfRansac(double errX, double errY, int nMinAcept, int nIntentos)
{
	L_Array<bool> compat;
	std::vector<double> xR, yR, xP, yP;
	double xRc[4], yRc[4], xPc[4], yPc[4];
	int i, j, nacept, azar;
	L_CalceNodo *cal, **p;
	L_TransfProyectivaPI2D invTr;
	double den, err;

	if ((int)calL.size() < nMinAcept)
		return false;

	if (ransac_2cal == true)
	{
		xR.resize(calL.size()*2);
		yR.resize(calL.size()*2);
		xP.resize(calL.size()*2);
		yP.resize(calL.size()*2);
	}
	else
	{
		xR.resize(calL.size());
		yR.resize(calL.size());
		xP.resize(calL.size());
		yP.resize(calL.size());
	}
	compat.resize((int)calL.size());

	for (j=0, cal=calL.root; j<(int)calL.size() && cal!=NULL; j++, cal=cal->sig)
	{
		if (distorsionarCoordImagenPrueba)
			camInter.distorsionaPixel(cal->c.dPru->x0, cal->c.dPru->y0, xP[j], yP[j]);
		else
			{xP[j] = cal->c.dPru->x0; yP[j] = cal->c.dPru->y0;}
		if (distorsionarCoordImagenReferencia)
			camInter.distorsionaPixel(cal->c.dRef->x0, cal->c.dRef->y0, xR[j], yR[j]);
		else
			{xR[j] = cal->c.dRef->x0; yR[j] = cal->c.dRef->y0;}
		if (ransac_2cal)
		{
			if (distorsionarCoordImagenPrueba)
				camInter.distorsionaPixel(cal->c.dPru->x0PuntaFlecha(), cal->c.dPru->y0PuntaFlecha(), xP[j+calL.size()], yP[j+calL.size()]);
			else
				{xP[j+calL.size()] = cal->c.dPru->x0PuntaFlecha(); yP[j+calL.size()] = cal->c.dPru->y0PuntaFlecha();}
			if (distorsionarCoordImagenReferencia)
				camInter.distorsionaPixel(cal->c.dRef->x0PuntaFlecha(), cal->c.dRef->y0PuntaFlecha(), xR[j+calL.size()], yR[j+calL.size()]);
			else
				{xR[j+calL.size()] = cal->c.dRef->x0PuntaFlecha(); yR[j+calL.size()] = cal->c.dRef->y0PuntaFlecha();}
		}
	}

	for (i=0; i<nIntentos; i++)
	{
		for (j=0; j<(int)calL.size(); j++)
			compat[j]=false;
		nacept=0;
		// Seleccionar calces iniciales
		nacept=0;
		while (nacept < 4)
		{
			azar=(int)( rand()%calL.size() );
			if (compat[azar]==true) // mala suerte
				continue;
			compat[azar]=true;
			xRc[nacept] = xR[azar];
			yRc[nacept] = yR[azar];
			xPc[nacept] = xP[azar];
			yPc[nacept] = yP[azar];
			nacept++;
			if (ransac_2cal)
			{
				xRc[nacept] = xR[azar+calL.size()];
				yRc[nacept] = yR[azar+calL.size()];
				xPc[nacept] = xP[azar+calL.size()];
				yPc[nacept] = yP[azar+calL.size()];
				nacept++;
			}
		}
		if (invTr.calcTransfExacto(xPc, yPc, xRc, yRc)==false) // Se llama con los argumentos al reves
			continue;
		for (j=0; j<(int)calL.size(); j++)
		{
			if (compat[j]==true)
				continue;
			// errores de retroproyeccion
			den = invTr.m31*xP[j] + invTr.m32*yP[j] + 1; //m33 = 1
			err=xR[j] - (invTr.m11*xP[j] + invTr.m12*yP[j] + invTr.tx)/den;
			if (err<0)
			{
				if (-err > errX)
					continue;
			}
			else
			{
				if (err > errX)
					continue;
			}
			err=yR[j] - (invTr.m21*xP[j] + invTr.m22*yP[j] + invTr.ty)/den;
			if (err<0)
			{
				if (-err > errY)
					continue;
			}
			else
			{
				if (err > errY)
					continue;
			}
			compat[j]=true;
			nacept++;
		}
		if (nacept>=nMinAcept)
		{
			p=&calL.root;
			for (j=0; j<(int)calL.size(); j++)
			{
				if (compat[j]==true)
					p=&(*p)->sig;
				else
				{
					cal = *p;
					*p = (*p)->sig;
					delete cal;
				}
			}
			calL.pult=p;
			*calL.pult=NULL;
			calL.n=nacept;
			calcTransfMinCuad();
			return true;
		}
	}
	return false;
}

bool L_TransfProyectivaPI2D::deformaImagen(L_ImageGrayDouble &dest, const L_ImageGrayDouble &orig, L_Transf2DDeformTipo tipo) const
{
	L_PuntInt p1, p2;
	L_TransfProyectivaPI2D itr;
	double x=0, y=0;
	int i;
	double punBorde[4][2]={ {0,0} , {1,0}, {1,1}, {0,1} };
	double nPunBorde=4;
	double coef[8];

	switch(tipo)
	{
	case PIXELADO:
	case BILIN_1NIV_ISOTR:
	case BILIN_XNIV_ISOTR:
		break;
	case BILIN_1NIV_ANISOTR:
	case BILIN_XNIV_ANISOTR:
	default:
		return false;
	}
	for (i=0; i<nPunBorde; i++)
	{
		p1.x0 = punBorde[i][0]*orig.lx;
		p1.y0 = punBorde[i][1]*orig.ly;
		proyeccionDe(p1, p2);
		if (p2.x0>x)
			x=p2.x0;
		if (p2.y0>y)
			y=p2.y0;
	}
	dest.reallocate((int)x, (int)y);
	if (!itr.calcParametrosInversosDe(*this))
		return false;
	switch(tipo)
	{
	case PIXELADO:
		for (i=0; i<8; i++)
			coef[i] = el(i);
		dest.AplicaTransfProyectivaPixelado(orig, coef);
		break;

	default:
		return false; // Hay que implementarlo
	}
	return true;
}

L_TransfProyectivaPI2D* L_TransfProyectivaPI2D::clona()
{
	int i;
	L_TransfProyectivaPI2D *ret = new L_TransfProyectivaPI2D();
	for (i=0; i<8; i++)
		ret->el(i) = el(i);
	ret->L_Transf2D_POD::operator=(*this);
	ret->calL = calL;
	return ret;
}

L_TransfProyectivaPI2D* L_TransfProyectivaPI2D::clonaRegalaCalces()
{
	int i;
	L_TransfProyectivaPI2D *ret = new L_TransfProyectivaPI2D();
	
	for (i=0; i<8; i++)
		ret->el(i) = el(i);
	ret->L_Transf2D_POD::operator=(*this);
	calL.moveListTo(ret->calL);
	return ret;
}




bool L_SURFGaussCalculador::calcDescriptor(L_ImageGrayDouble &grad1, L_ImageGrayDouble &grad2, bool gradPolares)
{
	double cosTheta, senTheta;
	double u, v; // posición relativa a la ventana de analisis
	int x, y; // posición absoluta en la imagen
	int nxt, nyt;
	L_LinearQuantizer c_i, c_j;

	long i, j; // indices
	int cont;

	double factorUV;
	double grmod, grang, grder, grarr;
	double angLocal;
	double cosAngThis, senAngThis;
	double xc, yc;
	int xIni, xFin, yIni, yFin, xrnd, yrnd;

	nt=nx*ny*4;

	//matriz de rotacion para "eje v" invertido
	cosTheta=cos(this->ang);
	senTheta=sin(this->ang);
	xc=this->pir_x-0.5; // Para pasar de this->x en la imagen normal a xc en el gradiente hacia adelante
	yc=this->pir_y-0.5; // Para pasar de this->x en la imagen normal a xc en el gradiente hacia adelante

	//definiciones de uso interno: para rotacion
	#if defined(calcU) || defined(calcV) || defined(u_) || defined(v_)
		#error #define conflicts
	#endif
	#define calcU(X,Y) L_ROUND( ((X)-xc)*cosTheta - ((Y)-yc)*senTheta ) // +0.5 para que quede redondeado
	#define calcV(X,Y) L_ROUND( ((X)-xc)*senTheta + ((Y)-yc)*cosTheta )
	#define u_ (L_FLOOR(u))
	#define v_ (L_FLOOR(v))

	nxt=lx*nx/2; // dg_nx DEBE ser par.
	nyt=ly*ny/2; // dg_ny DEBE ser par.

	xIni=L_FLOOR(xc-1.4142*nxt + 0.5);
	if (xIni<0)
		return false;
	xFin=L_FLOOR(xc+1.4142*nxt + 0.5);
	if (xFin>=grad1.lx)
		return false;
	yIni=L_FLOOR(yc-1.4142*nyt + 0.5);
	if (yIni<0)
		return false;
	yFin=L_FLOOR(yc+1.4142*nyt + 0.5);
	if (yFin>=grad1.ly)
		return false;

	for (cont=0; cont<nt; cont++)
		vector[cont]=0;
	if (gdx_h.size() == 0)
		L_calcGausHNorm(gdx_h,nxt,2*nxt); //gdx_h[i]: i debe ir entre (-nxt)+nxt y (nxt-1)+nxt
	if (gdy_h.size() == 0)
		L_calcGausHNorm(gdy_h,nyt,2*nyt); //gdx_h[i]: i debe ir entre (-nyt)+nyt y (nyt-1)+nyt

	c_i.setFromCellRange(-nxt,-nxt+lx,0);
	c_j.setFromCellRange(-nyt,-nyt+ly,0);

	cosAngThis = cos(-this->ang);
	senAngThis = sin(-this->ang);


	if (gradPolares)
	{
		for (y=yIni; y<=yFin; y++)
		{
			for (x=xIni; x<=xFin; x++)
			{
				u=calcU(x,y);
				v=calcV(x,y);

				if (u<-nxt || u>=nxt || v<-nyt || v>=nyt)
					continue;

				xrnd = L_ROUND(x);
				yrnd = L_ROUND(y);

				// Calculo del gradiente grmod<grang
				grmod = grad1.pix(xrnd,yrnd);
				grang = grad2.pix(xrnd,yrnd);
				angLocal=grang-this->ang;
				grder = grmod*cos(angLocal);
				grarr = grmod*sin(angLocal);
				// Calculo de i,j
				i=L_LinearQuantizer_nCluster(c_i,u);
				j=L_LinearQuantizer_nCluster(c_j,v);
				factorUV=gdx_h[u_+nxt]*gdy_h[v_+nyt];

				vector[i+j*nx + 0*nx*ny]+=factorUV*grder;
				vector[i+j*nx + 1*nx*ny]+=factorUV*grarr;
				vector[i+j*nx + 2*nx*ny]+=factorUV*L_ABS(grder);
				vector[i+j*nx + 3*nx*ny]+=factorUV*L_ABS(grarr);
			}
		}
	}
	else
	{
		for (y=yIni; y<=yFin; y++)
		{
			for (x=xIni; x<=xFin; x++)
			{
				u=calcU(x,y);
				v=calcV(x,y);

				if (u<-nxt || u>=nxt || v<-nyt || v>=nyt)
					continue;

				xrnd = L_ROUND(x);
				yrnd = L_ROUND(y);

				// Calculo del gradiente grmod<grang
				double grx = grad1.pix(xrnd,yrnd);
				double gry = grad2.pix(xrnd,yrnd);
				// z1 = grx + i*gry;  z2 = cosAngThis + i*senAngThis
				// z1*z2 = (grx*cosAngThis - gry*senAngThis) + i*(grx*senAngThis + gry*cosAngThis);
				grder = grx*cosAngThis - gry*senAngThis;
				grarr = grx*senAngThis + gry*cosAngThis;
				// Calculo de i,j
				i=L_LinearQuantizer_nCluster(c_i,u);
				j=L_LinearQuantizer_nCluster(c_j,v);
				factorUV=gdx_h[u_+nxt]*gdy_h[v_+nyt];

				vector[i+j*nx + 0*nx*ny]+=factorUV*grder;
				vector[i+j*nx + 1*nx*ny]+=factorUV*grarr;
				vector[i+j*nx + 2*nx*ny]+=factorUV*L_ABS(grder);
				vector[i+j*nx + 3*nx*ny]+=factorUV*L_ABS(grarr);
			}
		}
	}

	#undef calcU
	#undef calcV
	#undef u_
	#undef v_

	normalizaHist();
	for (cont=0; cont<nt; cont++)
	{
		if (vector[cont]>0.2)
			vector[cont]=0.2;
		if (vector[cont]<-0.2)
			vector[cont]=-0.2;
	}
	normalizaHist();
	return true;
}
