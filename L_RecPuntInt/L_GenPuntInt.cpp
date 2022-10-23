#include "L_GenPuntInt.h"

// Este archivo se puede compilar de modo optimizado aunque el codeMapping se compile en modo debug
// para hacer debug mas rapido

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


void L_GenDescrPuntInt::dibFlechasDescriptores(L_ShapeArray &lins)
{
	desL.dibujaFlechas(lins);
}

void L_GenDescrPuntInt::dibCuadradosDescriptores(L_ShapeArray &lins)
{
	desL.dibujaCuadrados(lins);
}

L_SSS::L_SSS():paramSSS("SSS",0)
{
	nio=2;
	nie=5;
	nsub=-1;
	imorig=NULL;
	s0=0.5;
	s1=1.6;
	s1_ND=1.0;
	sMin=1.2;
	reltam=L_Frac(0,1);
	dupl=false;
	submEfic=true;
	usasqrt2 = false;
	paramSSS.addFrom("nio",&nio);
	paramSSS.addFrom("nie",&nie);
	paramSSS.addFrom("s1",&s1);
	paramSSS.addFrom("s1_ND",&s1_ND);
	paramSSS.addFrom("sMin",&sMin);
	paramSSS.addFrom("dupl",&dupl);
	paramSSS.addFrom("submEfic", &submEfic);
	paramSSS.addFrom("usasqrt2", &usasqrt2);
}

L_SSS::L_SSS(int nio_, int nie_) : paramSSS("SSS",0)
{
	nio=nio_;
	nie=nie_;
	nsub=-1;
	imorig=NULL;
	s0=0.5;
	s1=1.6;
	s1_ND=1.0;
	sMin=1.2;
	reltam=L_Frac(0,1);
	dupl=false;
	submEfic=true;
	usasqrt2 = true;
	paramSSS.addFrom("nio",&nio);
	paramSSS.addFrom("nie",&nie);
	paramSSS.addFrom("s1",&s1);
	paramSSS.addFrom("s1_ND",&s1_ND);
	paramSSS.addFrom("sMin",&sMin);
	paramSSS.addFrom("dupl",&dupl);
	paramSSS.addFrom("submEfic", &submEfic);
	paramSSS.addFrom("usasqrt2", &usasqrt2);
}


void L_SSS::fijaImagenOriginal(L_ImageGrayDouble &im, double sigma)
{
	imorig=&im;
	s0=sigma;
}

void L_SSS::borraOctava()
{
	int i;
	for (i=0; i<nie; i++)
		imarr[i].destroy();
}

bool L_SSS::nuevaOctava()
{
	int iant, ides=-1;
	bool primVez=false;
	double sPrimImag;

	if (s0==-1)
		return false;
#ifdef L_BASIC_DEBUG
	if (nie<=0)
		L_hard_shutdown("L_SSS::nuevaOctava : non-positive size allocation");
#endif
	if (imarr.size() == 0)
	{
		imarr.resize(nie);
	}
	if (dupl)
		sPrimImag=s1;
	else
		sPrimImag=s1_ND;

	for (iant=nio; iant<nie; iant++)
	{
		ides=iant-nio;
		if (imarr[iant].data()!=NULL && (iant==nio || submEfic==true)) // submuestrear imarr[iant] a imarr[ides]
		{
			imarr[ides].subm(imarr[iant]);
		}
		else
		{
			if (iant==nio) // generar a partir de imorig la primera imagen del scale space.
			{
				primVez=true;
				nsub=0;
				if (sigmrel.size() == 0)
				{
					sigmrel.resize(nie);
					for (int i=0; i<nie; i++)
						sigmrel[i]=pow(2.0,i/(double)nio);
				}
				if (sigm.size() == 0)
				{
					sigm.resize(nie);
					for (int i=0; i<nie; i++)
						sigm[i]=sPrimImag*sigmrel[i];
				}
				if (dupl)
				{
					reltam=L_Frac(2,1);
					if (2*s0>=sigm[0])
						imarr[ides].doubleResolutionBilineal(*imorig);
					else
					{
						L_ImageGrayDouble imt;
						imt.doubleResolutionBilineal(*imorig);
						if (usasqrt2)
							imarr[ides].smoothXYSqrt2(imt);
						else
							imarr[ides].smooth_changeSigma(imt, 2*s0, sigm[0], sMin);
					}
				}
				else
				{
					reltam=L_Frac(1,1);
					if (s0>=sigm[0])
						imarr[ides]=*imorig;
					else
					{
						if (usasqrt2)
							imarr[ides].smoothXYSqrt2(*imorig);
						else
							imarr[ides].smooth_changeSigma(*imorig, s0, sigm[0], sMin);
					}
				}
			}
			else // generar a partir de imarr[ides-1]
			{
				if (usasqrt2)
					imarr[ides].smoothXYSqrt2(imarr[ides-1]);
				else
					imarr[ides].smooth_changeSigma(imarr[ides-1], sigm[ides-1], sigm[ides], sMin);
			}
		}
	}
	for (; ides<nie; ides++)
	{
		if (usasqrt2)
			imarr[ides].smoothXYSqrt2(imarr[ides-1]);
		else
			imarr[ides].smooth_changeSigma(imarr[ides-1], sigm[ides-1], sigm[ides], sMin);
	}
	if (!primVez)
	{
		L_Frac unmedio(1,2);
		reltam*=unmedio;
		nsub++;
	}
	return true;
};

void L_SSS::grabaArchivosBMP(const L_FileName &nomarch)
{
	L_FileName nom2;
	int num;
	int i;

	if (s0==-1)
		return;

	num=nsub*nie;

	for (i=0; i<nie; i++)
	{
		nom2=nomarch;
		nom2.name=nom2.name+L_String(num+i,2);
		imarr[i].writeBMP(nom2.c_str());
	}
}

void L_SSS::destroy()
{
}




L_PSS::L_PSS():paramPSS("PSS",0)
{
	porcSub=L_Frac(2,3);
	nsub=-1;
	s0=0.5;
	s1=1.4;
	s1_ND=1.0;
	sMin=1.2;
	reltam=L_Frac(0,1);
	dupl=false;
	usasqrt2=false;
	imprSigmas=false;

	paramPSS.addFrom("porcSub_1_num",&porcSub.num);
	paramPSS.addFrom("porcSub_2_den",&porcSub.den);
	paramPSS.addFrom("s1",&s1);
	paramPSS.addFrom("s1_ND",&s1_ND);
	paramPSS.addFrom("sMin",&sMin);
	paramPSS.addFrom("dupl",&dupl);
	paramPSS.addFrom("usasqrt2",&usasqrt2);
	paramPSS.addFrom("imprSigmas",&imprSigmas);
}
L_PSS::L_PSS(L_Frac &f):paramPSS("PSS",0)
{
	porcSub=f;
	nsub=-1;
	s0=0.5;
	s1=1.4;
	s1_ND=1.0;
	sMin=1.2;
	reltam=L_Frac(0,1);
	dupl=false;
	usasqrt2=false;
	imprSigmas=false;

	paramPSS.addFrom("reltam_1_num",&reltam.num);
	paramPSS.addFrom("reltam_2_den",&reltam.den);
	paramPSS.addFrom("s1",&s1);
	paramPSS.addFrom("s1_ND",&s1_ND);
	paramPSS.addFrom("sMin",&sMin);
	paramPSS.addFrom("dupl",&dupl);
	paramPSS.addFrom("usasqrt2",&usasqrt2);
	paramPSS.addFrom("imprSigmas",&imprSigmas);
}


void L_PSS::fijaImagenOriginal(L_ImageGrayDouble &im, double sigma)
{
	double sPrimImag;
	imorig=&im;
	s0=sigma;
	if (dupl)
		sPrimImag=s1;
	else
		sPrimImag=s1_ND;
	nsub=-1;
	if (dupl)
	{
		reltam=L_Frac(2,1);
		imantG.doubleResolutionBilineal(im);
		sigma*=2;
	}
	else
	{
		reltam=L_Frac(1,1);
		imantG=*imorig;
	}

	if (sigma<=sPrimImag)
	{
		if (usasqrt2)
		{
			imact.smoothXYSqrt2(imantG);
			if (imprSigmas)
				printf("LPSS::fijaImagenOriginal  sigma=sqrt(2)\n");
		}
		else
		{
			imact.smooth_changeSigma(imantG, sigma, sPrimImag, sMin);
			if (imprSigmas)
				printf("LPSS::fijaImagenOriginal  sigma=%.4f\n",sqrt(sPrimImag*sPrimImag-sigma*sigma));
		}
	}
	else
		imantG.swap(imact);
}


bool L_PSS::nuevaImagen()
{
	double sPrimImag;
	if (s0==-1)
		return false;
	if (dupl)
		sPrimImag=s1;
	else
		sPrimImag=s1_ND;
	// imact --> imantG -(smooth)-> imactG -(subm)-> imact
	imact.swap(imantG);
	if (usasqrt2)
	{
		imactG.smoothXYSqrt2(imantG);
		if (imprSigmas)
			printf("LPSS::nuevaImagen  sigma=sqrt(2)\n");
	}
	else
	{
		imactG.smooth_changeSigma(imantG, sPrimImag, sPrimImag/(double)porcSub);
		if (imprSigmas)
			printf("LPSS::fijaImagenOriginal  sigma=%.4f\n",sqrt(sPrimImag/(double)porcSub*sPrimImag/(double)porcSub-sPrimImag*sPrimImag));
	}
	imact.subm(imactG, porcSub);
	if (nsub!=-1)
		reltam*=porcSub;
	nsub++;
	return true;
}

void L_PSS::grabaArchivosBMP(const L_FileName &nomarch)
{
	L_FileName nom2;
	nom2=nomarch;
	nom2.name=nom2.name+L_String(nsub,2);
	imact.writeBMP(nom2.c_str());
}

void L_STrioEstat::calc_M_BTB_BT(L_Matrix &M_BTB_BT)
{
	L_Matrix M_B;
	L_Matrix M_BT;
	L_Matrix M_BTB;
	int x,y,n;
	int ind=0;
	M_B.reallocate(27,10);
	for (n=-1; n<=1; n++)
	{
		for (x=-1; x<=1; x++)
		{
			for (y=-1; y<=1; y++)
			{
				M_B(ind,0)=x*x;
				M_B(ind,1)=y*y;
				M_B(ind,2)=n*n;
				M_B(ind,3)=x*y;
				M_B(ind,4)=x*n;
				M_B(ind,5)=y*n;
				M_B(ind,6)=x;
				M_B(ind,7)=y;
				M_B(ind,8)=n;
				M_B(ind,9)=1;
				ind++;
			}
		}
	}
	M_BT.transpOf(M_B);
	M_BTB.OP_mult(M_BT,M_B);
	M_BTB.invertMe();
	M_BTB_BT.OP_mult(M_BTB,M_BT);
}

void L_PTrioEstat::calc_M_BTB_BT(L_Matrix &M_BTB_BT)
{
	L_Matrix M_B;
	L_Matrix M_BT;
	L_Matrix M_BTB;
	int x,y;
	int ind=0;
	M_B.reallocate(11,8);
	x=0; y=0;
	M_B(ind,0)=0;  //x*x
	M_B(ind,1)=0;  //y*y
	M_B(ind,2)=1;  //n*n;
	M_B(ind,3)=0;  //x*y
	M_B(ind,4)=0;  //x
	M_B(ind,5)=0;  //y
	M_B(ind,6)=-1; //n
	M_B(ind,7)=1;
	ind++;
	for (x=-1; x<=1; x++)
	{
		for (y=-1; y<=1; y++)
		{
			M_B(ind,0)=x*x;
			M_B(ind,1)=y*y;
			M_B(ind,2)=0;
			M_B(ind,3)=x*y;
			M_B(ind,4)=x;
			M_B(ind,5)=y;
			M_B(ind,6)=0;
			M_B(ind,7)=1;
			ind++;
		}
	}
	M_B(ind,0)=0;  //x*x
	M_B(ind,1)=0;  //y*y
	M_B(ind,2)=1;  //n*n;
	M_B(ind,3)=0;  //x*y
	M_B(ind,4)=0;  //x
	M_B(ind,5)=0;  //y
	M_B(ind,6)=1; //n
	M_B(ind,7)=1;
	ind++;
	M_BT.transpOf(M_B);
	M_BTB.OP_mult(M_BT,M_B);
	M_BTB.invertMe();
	M_BTB_BT.OP_mult(M_BTB,M_BT);
}

void L_STrioIm::ajustaParabol(int i, int j, L_Matrix &grad, L_Matrix &hess)
{
	int x, y, n, ind=0;
	L_Matrix fn;
	L_Matrix a;
	fn.reallocate(27,1);

	for (n=-1; n<=1; n++)
	{
		for (y=-1; y<=1; y++)
		{
			for (x=-1; x<=1; x++)
			{
				fn(ind,0)=im[n+1].pix(i+x,j+y);
				ind++;
			}
		}
	}
	if (M_BTB_BT.begin() == NULL)
		L_STrioEstat::calc_M_BTB_BT(M_BTB_BT);
	L_Matrix M_fn;
	M_fn.OP_mult(M_BTB_BT , fn);
	M_fn.swap(a);
	hess.reallocate(3,3);
	grad.reallocate(3,1);
	hess(0,0)=a(0,0)*2;// coef (x*x)
	hess(1,0)=a(3,0);  // coef (x*y)
	hess(2,0)=a(4,0);  // coef (x*n)
	hess(0,1)=a(3,0);  // coef (x*y)
	hess(1,1)=a(1,0)*2;// coef (y*y)
	hess(2,1)=a(5,0);  // coef (y*n)
	hess(0,2)=a(4,0);  // coef (x*n)
	hess(1,2)=a(5,0);  // coef (y*n)
	hess(2,2)=a(2,0)*2;// coef (n*n)
	grad(0,0)=a(6,0);  // coef (x)
	grad(1,0)=a(7,0);  // coef (y)
	grad(2,0)=a(8,0);  // coef (n)
}

void L_STrioIm::ajustaTaylor1(int i, int j, L_Matrix &grad, L_Matrix &hess)
{
	hess.reallocate(3,3);
	grad.reallocate(3,1);

	grad(0,0)=(im[2].pix(i,j)-im[0].pix(i,j))/2;
	grad(1,0)=(im[1].pix(i+1,j)-im[1].pix(i-1,j))/2;
	grad(2,0)=(im[1].pix(i,j+1)-im[1].pix(i,j-1))/2;

	hess(0,0)=im[2].pix(i,j)+im[0].pix(i,j)-2*im[1].pix(i,j);
	hess(1,1)=im[1].pix(i+1,j)+im[1].pix(i-1,j)-2*im[1].pix(i,j);
	hess(2,2)=im[1].pix(i,j+1)+im[1].pix(i,j-1)-2*im[1].pix(i,j);

	hess(0,1)=hess(1,0)=(
		(im[2].pix(i+1,j)-im[0].pix(i+1,j))-
		(im[2].pix(i-1,j)-im[0].pix(i-1,j))
	)/4;

	hess(0,2)=hess(2,0)=(
		(im[2].pix(i,j+1)-im[0].pix(i,j+1))-
		(im[2].pix(i,j-1)-im[0].pix(i,j-1))
	)/4;

	hess(1,2)=hess(2,1)=(
		(im[1].pix(i+1,j+1)-im[1].pix(i-1,j+1))-
		(im[1].pix(i+1,j-1)-im[1].pix(i-1,j-1))
	)/4;
}

double L_STrioIm::interpola(L_PuntInt &p,  L_Matrix &xmax, const L_Matrix &grad, const L_Matrix &gradtr, const L_Matrix& hess, const L_Matrix &hessinv)
{
	double val;
	L_Matrix pr;
	// "(-hessinv*grad).swap(xmax);"
	{
		pr.OP_mult(hessinv,grad);
		xmax.OP_mult(pr, -1);
	}
	val=im[1].pix((int)p.pir_x,(int)p.pir_y);

    p.pir_y+=xmax(1,0);

	if (xmax(0,0)<-0.5)
	{
		if (xmax(0,0)<-1)
			xmax(0,0)=-0.99;
		p.pir_x-=xmax(0,0);
	}
	else if (xmax(0,0)>0.5)
	{
		if (xmax(0,0)>1)
			xmax(0,0)=0.99;
		p.pir_x+=xmax(0,0);
	}
	if (xmax(1,0)<-0.5)
	{
		if (xmax(1,0)<-1)
			xmax(1,0)=-0.99;
		p.pir_y-=xmax(1,0);
	}
	else if (xmax(1,0)>0.5)
	{
		if (xmax(1,0)>1)
			xmax(1,0)=0.99;
		p.pir_y+=xmax(1,0);
	}
	if (xmax(2,0)<-0.5)
	{
		if (xmax(2,0)<1)
			xmax(2,0)=-0.99;
		p.pir_n--;
	}
	else if (xmax(2,0)>0.5)
	{
		if (xmax(2,0)>1)
			xmax(2,0)=0.99;
		p.pir_n++;
	}
	L_Matrix gm;
	gm.OP_mult(gradtr,xmax);
	val=val+0.5*gm(0,0);
	return val;
}

void L_PTrioIm::ajustaParabol(int i, int j, L_Matrix &grad, L_Matrix &hess)
{
	int x, y, ind=0;
	L_Matrix fn;
	L_Matrix a;
	fn.reallocate(11,1);

	fn(ind,0)=im[0]->pix(i*f01->den/f01->num,j*f01->den/f01->num);
	ind++;

	for (y=-1; y<=1; y++)
	{
		for (x=-1; x<=1; x++)
		{
			fn(ind,0)=im[1]->pix(i+x,j+y);
			ind++;
		}
	}

	x=0; y=0; n=1;
	fn(ind,0)=im[0]->pix(i*f01->num/f01->den,j*f01->num/f01->den);
	ind++;

	if (M_BTB_BT.begin() == NULL)
		L_STrioEstat::calc_M_BTB_BT(M_BTB_BT);
	L_Matrix Mfn;
	Mfn.OP_mult(M_BTB_BT , fn);
	Mfn.swap(a);
	hess.reallocate(3,3);
	grad.reallocate(3,1);
	hess(0,0)=a(0,0)*2;// coef (x*x)
	hess(1,0)=a(3,0);  // coef (x*y)
	hess(2,0)=0;             // coef (x*n)
	hess(0,1)=a(3,0);  // coef (x*y)
	hess(1,1)=a(1,0)*2;// coef (y*y)
	hess(2,1)=0;             // coef (y*n)
	hess(0,2)=0;             // coef (x*n)
	hess(1,2)=0;             // coef (y*n)
	hess(2,2)=a(2,0)*2;// coef (n*n)
	grad(0,0)=a(4,0);  // coef (x)
	grad(1,0)=a(5,0);  // coef (y)
	grad(2,0)=a(6,0);  // coef (n)
}

void L_PTrioIm::ajustaTaylor1(int i, int j, L_Matrix &grad, L_Matrix &hess)
{
	hess.reallocate(3,3);
	grad.reallocate(3,1);

	grad(0,0)=(im[2]->pix(i,j)-im[0]->pix(i,j))/2;
	grad(1,0)=(im[1]->pix(i+1,j)-im[1]->pix(i-1,j))/2;
	grad(2,0)=(im[1]->pix(i,j+1)-im[1]->pix(i,j-1))/2;

	hess(0,0)=im[2]->pix(i,j)+im[0]->pix(i,j)-2*im[1]->pix(i,j);
	hess(1,1)=im[1]->pix(i+1,j)+im[1]->pix(i-1,j)-2*im[1]->pix(i,j);
	hess(2,2)=im[1]->pix(i,j+1)+im[1]->pix(i,j-1)-2*im[1]->pix(i,j);

	hess(0,1)=hess(1,0)=0; // No se puede estimar bien en la piramide
	hess(0,2)=hess(2,0)=0; // No se puede estimar bien en la piramide
	hess(1,2)=hess(2,1)=0; // No se puede estimar bien en la piramide
}


double L_PTrioIm::interpola(L_PuntInt &p,  L_Matrix &xmax, const L_Matrix &grad, const L_Matrix &gradtr, const L_Matrix& hess, const L_Matrix &hessinv)
{
	double val;
	val=im[1]->pix((int)p.pir_x,(int)p.pir_y);
	// "(-hessinv*grad).swap(xmax);"
	L_Matrix rt;
	rt.OP_mult(hessinv , grad);
	rt.OP_mult(-1);
	rt.swap(xmax);

	p.pir_x+=xmax(0,0);
    p.pir_y+=xmax(1,0);

	if (xmax(0,0)<-0.5)
	{
		if (xmax(0,0)<-1)
			return 0.0; // botar el punto
	}
	else if (xmax(0,0)>0.5)
	{
		if (xmax(0,0)>1)
			return 0.0; // botar el punto
	}
	if (xmax(1,0)<-0.5)
	{
		if (xmax(1,0)<-1)
			return 0.0; // botar el punto
	}
	else if (xmax(1,0)>0.5)
	{
		if (xmax(1,0)>1)
			return 0.0; // botar el punto
	}
	if (xmax(2,0)<-0.5)
    {
			return 0.0; // botar el punto
    }
	else if (xmax(2,0)>0.5)
    {
		    return 0.0; // botar el punto
    }
	//"val=val+0.5*(gradtr*xmax)(0,0);"
	{
		L_Matrix gx;
		gx.OP_mult(gradtr, xmax);
		val=val+0.5*gx(0,0);
	}
	return val;
}

////////////////////////////
//////// DoG submuestreado
///////////////////////////

L_SDoG_SIFT::L_SDoG_SIFT() : paramSDOGSIFT("SDOGSIFT",5)
{
	tipoP = SDOG;
	tipoD = SIFT;
	r_lin=10;
	val_min_DoG=0.001; //0.01;
	usa_parabol=2; nobj=-1;
	minX=30; minY=30; // de L_GenDescrPuntInt
	porcCurvMax=0.94;
	s0_def=0.3; s1_def=1.26;  s1_def_chico=1.26;
	calces_maxDifIntens=0.06;
	ultId = 0;
	// Parametros de la class_name padre
	paramSDOGSIFT.addFrom("minX",&minX);
	paramSDOGSIFT.addFrom("minY",&minY);
	paramSDOGSIFT.addFrom("calces_maxDist",&calces_maxDist);
	paramSDOGSIFT.addFrom("calces_maxRatio",&calces_maxRatio);
	paramSDOGSIFT.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramSDOGSIFT.addFrom("calces_probCorrecto",&calces_probCorrecto);
	paramSDOGSIFT.addFrom("calces_maxDifIntens",&calces_maxDifIntens);
	// Parámetros de esta class_name
	paramSDOGSIFT.addFrom("r_lin",&r_lin);
	paramSDOGSIFT.addFrom("val_min_DoG",&val_min_DoG);
	paramSDOGSIFT.addFrom("val_minP_DoG",&val_minP_DoG);
	paramSDOGSIFT.addFrom("usa_parabol",&usa_parabol);
	paramSDOGSIFT.addFrom("porcCurvMax", &porcCurvMax);
	paramSDOGSIFT.addFrom("s0_def",&s0_def);
	paramSDOGSIFT.addFrom("s1_def",&s1_def);
	paramSDOGSIFT.addFrom("s1_def_chico",&s1_def_chico);

	paramSDOGSIFT.addChildren(&(paramSSS));
	paramSDOGSIFT.addChildren(&sgen.paramSIFT);
	paramSDOGSIFT.addChildren(&sdirL.paramSIFTDL);
}


bool L_SDoG_SIFT::calcPuntInt_act()
{
	int i, j;
	L_PuntInt pin;
	pin.hasDepth = 0;
	pin.tipoP = tipoP;

	double dx, dy, ds;
	pin.tipoP=SDOG;
	L_Matrix xmax, grad, gradtr;
	L_Matrix hess, hessXY, hessinv;
	double trazaXY,detXY,val;
	xmax.reallocate(3,1);
	xmax(0,0)=0; xmax(1,0)=0; xmax(2,0)=0;

	for (j=4; j<trioSDoG.im[1].ly-4; j++)
	{
		for (i=4; i<trioSDoG.im[1].lx-4; i++)
		{
			if (
				(trioSDoG.esMax3D(i,j) && trioSDoG.im[1].pix(i,j) >= val_min_DoG && trioSDoG.im[1].esMax2D_anillo(trioSDoG.im[1].pix(i,j)*porcCurvMax, i, j, 3*sigmrel[it]))
				||
				(trioSDoG.esMin3D(i,j) && trioSDoG.im[1].pix(i,j) <= -val_min_DoG && trioSDoG.im[1].esMin2D_anillo(trioSDoG.im[1].pix(i,j)*porcCurvMax, i, j, 3*sigmrel[it]))
				)
			{
				pin.pir_x=i;
				pin.pir_y=j;
				pin.pir_n=trioSDoG.n + 1; // trio 0 => imagen central 1
				pin.imorig=imorig;
				pin.imorig_lx=imorig->lx;
				pin.imorig_ly=imorig->ly;
				pin.hasDepth = 0;
				pin.pir_reltam=reltam;
				pin.intens=fabs(trioSDoG.im[1].pix(i,j));
				pin.signo = trioSDoG.im[1].pix(i,j) > 0 ? 1 : -1;

				switch (usa_parabol)
				{
				case 0:
					trioSDoG.ajustaTaylor1((int)pin.pir_x,(int)pin.pir_y,grad,hess);
					detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
					trazaXY=hess(0,0)+hess(1,1);
					gradtr.transpOf(grad);
					hessinv=hess; hessinv.invertMe();
					val=trioSDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
					break;
				case 1:
					trioSDoG.ajustaParabol((int)pin.pir_x,(int)pin.pir_y,grad,hess);
					detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
					trazaXY=hess(0,0)+hess(1,1);
					gradtr.transpOf(grad);
					hessinv=hess; hessinv.invertMe();
					val=trioSDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
					break; // Error : faltaba esto
				case 2:
					dx = L_INTERP_MAX_PARABOL(trioSDoG.im[1].pix(i-1,j), trioSDoG.im[1].pix(i,j), trioSDoG.im[1].pix(i+1,j));
					dy = L_INTERP_MAX_PARABOL(trioSDoG.im[1].pix(i,j-1), trioSDoG.im[1].pix(i,j), trioSDoG.im[1].pix(i,j+1));
					ds = L_INTERP_MAX_PARABOL(trioSDoG.im[0].pix(i,j), trioSDoG.im[1].pix(i,j), trioSDoG.im[2].pix(i,j));
					pin.pir_x += dx;
					pin.pir_y += dy;
					xmax(2,0) = ds;
					val = trioSDoG.im[1].pix(i,j);
					break;
				}

				pin.x0=pin.pir_x/(double)reltam; // Posicion relativa a la imagen original
				pin.y0=pin.pir_y/(double)reltam; // Posicion relativa a la imagen original
				
				pin.sigma0 = (dupl ? s1 : s1_ND)*pow(2.0,(trioSDoG.n+xmax(2,0))/(double)nio);

				punL.push_back(pin);
			}
		}
	}
	return true;
}

bool L_SDoG_SIFT::calcDirecc_act()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n-trioSDoG.n; // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, trioGradR.im[nim], trioGradA.im[nim], true, sigm[nim], ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_SDoG_SIFT::calcDescr_act()
{
	L_DireccNodo *p;
	int nim;
	int n=0;

	sgen.nobj=nobj;
	sgen.refGen=this;
	sgen.tipoP = tipoP;
	sgen.tipoD = tipoD;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n-trioSDoG.n;  // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sgen.L_Descriptor::operator=( p->c );
		sgen.id = ultId++;
		if (sgen.calcDescriptor(trioGradR.im[nim], trioGradA.im[nim], true))
			desL.push_back(sgen);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}

bool L_SDoG_SIFT::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	this->s1=s1_def;
	this->s1_ND=s1_def_chico;
	if (s0<0.1)
		s0=s0_def;
	fijaImagenOriginal(im, s0);
	nuevaOctava();
	while (1)
	{
		for (it=0; it<nio; it++)
		{
			if (it==0)
			{
				trioSDoG.im[0].OP_subtract(imarr[1],imarr[0]);
				trioSDoG.im[1].OP_subtract(imarr[2],imarr[1]);
				trioSDoG.im[2].OP_subtract(imarr[3],imarr[2]);
				trioGradR.im[0].gradR_2(trioGradA.im[0], imarr[0]);
				trioGradR.im[1].gradR_2(trioGradA.im[1], imarr[1]);
				trioGradR.im[2].gradR_2(trioGradA.im[2], imarr[2]);
			}
			else
			{
				trioSDoG.im[1].swap(trioSDoG.im[0]);
				trioSDoG.im[2].swap(trioSDoG.im[1]);
				trioSDoG.im[2].OP_subtract(imarr[it+3],imarr[it+2]);

				trioGradR.im[1].swap(trioGradR.im[0]);
				trioGradR.im[2].swap(trioGradR.im[1]);
				trioGradA.im[1].swap(trioGradA.im[0]);
				trioGradA.im[2].swap(trioGradA.im[1]);
				trioGradR.im[2].gradR_2(trioGradA.im[2], imarr[it+2]);
			}
			trioSDoG.n=nsub*nio+it;
			calcPuntInt_act();
			calcDirecc_act();
			calcDescr_act();
		}
		if (this->imarr[nie-1].lx<=minX || this->imarr[nie-1].ly<=minY)
			break;
		//if (this->imarr[nie-1].lx<=6*this->sigm[nie-1] || this->imarr[nie-1].ly<=6*this->sigm[nie-1]) // Imagen mas chica que el kernel
		//	break;
		nuevaOctava();
	}
	borraOctava();
	return true;
}

L_ParamManagerLocal *L_SDoG_SIFT::pideParams()
{
	return &paramSDOGSIFT;
}



////////////////////////////////////////////////////////////////////
// DoG piramide

L_PDoG_SIFT::L_PDoG_SIFT() : paramPDOGSIFT("PDOGSIFT",5)
{
	tipoP = PDOG;
	tipoD = SIFT;
	r_lin=10;
	val_minP_DoG=0.01/255;
	usa_parabol=2;
	nobj=-1;
	minX=30;
	minY=30; // de L_GenDescrPuntInt
	nPisos=0; // Pisos de la pirámide de gauss
	porcCurvMax=0.94;
	s0_def=0.3;
	s1_def=1.26;
	s1_def_chico=1.26;
	calces_maxDifIntens=0.06;
	ultId = 0;
	// Parametros de la class_name padre
	paramPDOGSIFT.addFrom("minX",&minX);
	paramPDOGSIFT.addFrom("minY",&minY);
	paramPDOGSIFT.addFrom("calces_maxDist",&calces_maxDist);
	paramPDOGSIFT.addFrom("calces_maxRatio",&calces_maxRatio);
	paramPDOGSIFT.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramPDOGSIFT.addFrom("calces_probCorrecto",&calces_probCorrecto);
	// Parametros de esta class_name
	paramPDOGSIFT.addFrom("r_lin",&r_lin);
	paramPDOGSIFT.addFrom("val_minP_DoG",&val_minP_DoG);
	paramPDOGSIFT.addFrom("usa_parabol",&usa_parabol);
	paramPDOGSIFT.addFrom("porcCurvMax", &porcCurvMax);
	paramPDOGSIFT.addFrom("s0_def",&s0_def);
	paramPDOGSIFT.addFrom("s1_def",&s1_def);
	paramPDOGSIFT.addFrom("s1_def_chico",&s1_def_chico);
	trioPDoG.f01=&genPiram.porcSub;
	trioPDoG.f12=&genPiram.porcSub;

	paramPDOGSIFT.addChildren(&genPiram.paramPSS);
	paramPDOGSIFT.addChildren(&sgen.paramSIFT);
	paramPDOGSIFT.addChildren(&sdirL.paramSIFTDL);
}

bool L_PDoG_SIFT::calcPiramides(L_ImageGrayDouble &im, double s0)
{
	int it=0;

	genPiram.fijaImagenOriginal(im, s0);
	while (genPiram.nuevaImagen() == true)
	{
		pirLaplace[it]=genPiram.imantG;
		if (pirLaplace[it].lx<10 || pirLaplace[it].ly<10)
			break;
		pirLaplace[it]-=genPiram.imactG;
		reltam[it]=genPiram.reltam;
		pirGauss[it]=genPiram.imactG;
		pirGradR[it].gradR_2(pirGradA[it], genPiram.imantG);

		it++;
		if (it>=L_PDS_NPIR)
			break;
	}
	nPisos=it;
	return true;
}

bool L_PDoG_SIFT::calcRegiones(L_ImageGrayDouble &im)
{
	int it, i, j;
	L_PuntInt pin;
	pin.imorig = &im;
	pin.imorig_lx = im.lx;
	pin.imorig_ly = im.ly;
	pin.hasDepth = 0;
	pin.tipoP = tipoP;

	L_Matrix xmax, grad, gradtr;
	L_Matrix hess, hessXY, hessinv;
	double trazaXY,detXY,val;
	double dx, dy, ds;
	
	xmax.reallocate(3,1);
	xmax(0,0)=0; xmax(1,0)=0; xmax(2,0)=0;


	for (it=1; it<nPisos-1; it++)
	{
		trioPDoG.im[0]=&pirLaplace[it-1];
		trioPDoG.im[1]=&pirLaplace[it];
		trioPDoG.im[2]=&pirLaplace[it+1];

		pin.pir_n=it;
		pin.sigma0=sigm[it-1];
		pin.pir_reltam=reltam[it];

		for (j=4; j<pirLaplace[it].ly-4; j++)
		{
			for (i=4; i<pirLaplace[it].lx-4; i++)
			{
				if (
					(trioPDoG.im[1]->pix(i,j) > val_minP_DoG && L_PTrioIm_esMax3D(trioPDoG,i,j) && trioPDoG.im[1]->esMax2D_anillo(trioPDoG.im[1]->pix(i,j)*porcCurvMax, i, j, 3))
					||
					(trioPDoG.im[1]->pix(i,j) < -val_minP_DoG && L_PTrioIm_esMin3D(trioPDoG,i,j) && trioPDoG.im[1]->esMin2D_anillo(trioPDoG.im[1]->pix(i,j)*porcCurvMax, i, j, 3))
					)
				{
					pin.pir_x=i;
					pin.pir_y=j;
					pin.x0=i/(double)reltam[it];
					pin.y0=j/(double)reltam[it];
					pin.intens=fabs(trioPDoG.im[1]->pix(i,j));
					pin.signo = trioPDoG.im[1]->pix(i,j) > 0 ? 1 : -1;

					switch (usa_parabol)
					{
					case 0:
						trioPDoG.ajustaTaylor1((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break;
					case 1:
						trioPDoG.ajustaParabol((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break; // Error : faltaba esto
					case 2:
						dx = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i-1,j), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i+1,j));
						dy = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i,j-1), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i,j+1));
						ds = L_INTERP_MAX_PARABOL(trioPDoG.pixTrio(0,i,j), trioPDoG.pixTrio(1,i,j), trioPDoG.pixTrio(2,i,j));
						pin.pir_x += dx;
						pin.pir_y += dy;
						xmax(2,0) = ds;
						val = trioPDoG.im[1]->pix(i,j);
						break;
					}

					pin.x0=pin.pir_x/(double)reltam[it];
					pin.y0=pin.pir_y/(double)reltam[it];
					pin.intens = fabs(val);
					pin.signo = val > 0 ? 1 : -1;

					pin.sigma0=(genPiram.dupl ? genPiram.s1 : genPiram.s1_ND) *pow((double)genPiram.porcSub,(double)(-(it+ds)));

					punL.push_back(pin);
				}
			}
		}
	}
	return true;
}

bool L_PDoG_SIFT::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n;
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, pirGradR[nim], pirGradA[nim], true, sigm[0], ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_PDoG_SIFT::calcDescriptores()
{
	L_DireccNodo *p;
	int nim;
	int n=0;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n;  // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sgen.L_Descriptor::operator=( p->c );

		sgen.nobj=nobj;
		sgen.refGen=this;
		sgen.tipoP = tipoP;
		sgen.tipoD = tipoD;

		sgen.id = ultId++;
		if (sgen.calcDescriptor(pirGradR[nim], pirGradA[nim], true))
			desL.push_back(sgen);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}

bool L_PDoG_SIFT::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	int i;
	genPiram.s1=s1_def;
	genPiram.s1_ND=s1_def_chico;
	if (s0<0.1)
		s0=s0_def;

	if (genPiram.dupl)
		for (i=0; i<L_PDS_NPIR; i++)
			sigm[i]=genPiram.s1*pow((double)genPiram.porcSub,(double)(-i));
	else
		for (i=0; i<L_PDS_NPIR; i++)
			sigm[i]=genPiram.s1_ND*pow((double)genPiram.porcSub,(double)(-i));
	calcPiramides(im, s0);
	calcRegiones(im);
	calcDirecciones();
	calcDescriptores();
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////
//// Harris submuestreado
/////////////////////////


L_SHarris_SIFT::L_SHarris_SIFT():paramSHARRISSIFT("SHARRISSIFT",3)
{
	tipoP = SHARRIS;
	tipoD = SIFT;
	factorS=1.5;
	alfa=0.04;
	val_minP_Harris=0.0001;
	val_minP_DoG=0.0003/2;
	nobj=-1;
	minX=30;
	minY=30;  // de L_GenDescrPuntInt
	porcCurvMax=0.94;
	s0_def=1.0;
	s1_def=1.6;
	s1_def_chico=1.4;
	usa_parabol = 2;
	calces_maxDifIntens=0.06;
	// Parametros de la class_name padre
	paramSHARRISSIFT.addFrom("minX",&minX);
	paramSHARRISSIFT.addFrom("minY",&minY);
	paramSHARRISSIFT.addFrom("calces_maxDist",&calces_maxDist);
	paramSHARRISSIFT.addFrom("calces_maxRatio",&calces_maxRatio);
	paramSHARRISSIFT.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramSHARRISSIFT.addFrom("calces_probCorrecto",&calces_probCorrecto);
	paramSHARRISSIFT.addFrom("calces_maxDifIntens",&calces_maxDifIntens);
	// Parametros de esta class_name
	paramSHARRISSIFT.addFrom("factorS",&factorS);
	paramSHARRISSIFT.addFrom("alfa",&alfa);
	paramSHARRISSIFT.addFrom("val_minP_DoG",&val_minP_DoG);
	paramSHARRISSIFT.addFrom("val_minP_Harris",&val_minP_Harris);
	paramSHARRISSIFT.addFrom("porcCurvMax", &porcCurvMax);
	paramSHARRISSIFT.addFrom("s0_def",&s0_def);
	paramSHARRISSIFT.addFrom("s1_def",&s1_def);
	paramSHARRISSIFT.addFrom("s1_def_chico",&s1_def_chico);
	paramSHARRISSIFT.addFrom("usa_parabol",&usa_parabol);

	paramSHARRISSIFT.addChildren(&(paramSSS));
	paramSHARRISSIFT.addChildren(&sgen.paramSIFT);
	paramSHARRISSIFT.addChildren(&sdirL.paramSIFTDL);
}

bool L_SHarris_SIFT::calcPuntInt_act()
{
	int i, j;
	L_PuntInt pin;
	pin.tipoP = SHARRIS;
	pin.hasDepth = 0;
	L_Matrix xmax, grad, gradtr;
	L_Matrix hess, hessXY, hessinv;
	double trazaXY,detXY,val;
	double dx, dy, ds;
	L_ImageGrayDouble imHarr;

	imHarr.harris(imarr[it+2], sigm[it+2], factorS*sigm[it+2], alfa);

	for (j=4; j<trioSDoG.im[1].ly-4; j++)
	{
		for (i=4; i<trioSDoG.im[1].lx-4; i++)
		{
			if
			(
				imHarr.esMax2D(i,j) && imHarr.pix(i,j)>=val_minP_Harris
				&&
				(
					( trioSDoG.esMax1D(i,j)&&trioSDoG.im[1].pix(i,j)>=val_minP_DoG )
					||
					( trioSDoG.esMin1D(i,j)&&trioSDoG.im[1].pix(i,j)<=-val_minP_DoG )
				)
				&&
				imHarr.esMax2D_anillo(imHarr.pix(i,j)*porcCurvMax, i, j, 5*sigmrel[it])
			)
			{
				pin.pir_x=i;
				pin.pir_y=j;
				pin.pir_n=it; // trio 0 => imagen central 1

				switch (usa_parabol)
				{
				case 0:
					trioSDoG.ajustaTaylor1((int)pin.pir_x,(int)pin.pir_y,grad,hess);
					detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
					trazaXY=hess(0,0)+hess(1,1);
					gradtr.transpOf(grad);
					hessinv=hess; hessinv.invertMe();
					ds = xmax(2,0);
					val=trioSDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
					break;
				case 1:
					trioSDoG.ajustaParabol((int)pin.pir_x,(int)pin.pir_y,grad,hess);
					detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
					trazaXY=hess(0,0)+hess(1,1);
					gradtr.transpOf(grad);
					hessinv=hess; hessinv.invertMe();
					ds = xmax(2,0);
					val=trioSDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
					break; // Error: faltaba esto
				case 2:
					dx = L_INTERP_MAX_PARABOL(trioSDoG.im[1].pix(i-1,j), trioSDoG.im[1].pix(i,j), trioSDoG.im[1].pix(i+1,j));
					dy = L_INTERP_MAX_PARABOL(trioSDoG.im[1].pix(i,j-1), trioSDoG.im[1].pix(i,j), trioSDoG.im[1].pix(i,j+1));
					ds = L_INTERP_MAX_PARABOL(trioSDoG.im[0].pix(i,j), trioSDoG.im[1].pix(i,j), trioSDoG.im[2].pix(i,j));
					pin.pir_x += dx;
					pin.pir_y += dy;
					xmax(2,0) = ds;
					val = trioSDoG.im[1].pix(i,j);
					break;
				}

				pin.intens=fabs(val);
				pin.signo = val > 0 ? 1 : -1;
				pin.imorig=imorig;
				pin.imorig_lx=imorig->lx;
				pin.imorig_ly=imorig->ly;
				pin.pir_reltam=reltam;

				pin.x0=pin.pir_x/(double)reltam;
				pin.y0=pin.pir_y/(double)reltam;
				pin.sigma0=(dupl ? s1 : s1_ND)*pow(2.0,(pin.pir_n+ds)/(double)nio);

				punL.push_back(pin);
			}
		}
	}
	return true;
}

bool L_SHarris_SIFT::calcDirecc_act()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n-trioSDoG.n;
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, trioGradR.im[nim], trioGradA.im[nim], true, sigm[nim], ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_SHarris_SIFT::calcDescr_act()
{
	L_DireccNodo *p;
	int nim;
	int n=0;
	sgen.nobj=nobj;
	sgen.refGen=this;
	sgen.tipoP = tipoP;
	sgen.tipoD = tipoD;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n-trioSDoG.n;
		sgen.L_Descriptor::operator=( p->c );
		sgen.id = ultId++;
		if (sgen.calcDescriptor(trioGradR.im[nim], trioGradA.im[nim], true))
		{
			desL.push_back(sgen);
		}
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}

bool L_SHarris_SIFT::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	this->s1=s1_def;
	this->s1_ND=s1_def_chico;
	if (s0<0.1)
		s0=s0_def;

	fijaImagenOriginal(im, s0);
	nuevaOctava();
	while (1)
	{
		for (it=0; it<nio; it++)
		{
			trioSDoG.n=nsub*nio+it;
			if (it==0)
			{
				trioSDoG.im[0].OP_subtract(imarr[1],imarr[0]);
				trioSDoG.im[1].OP_subtract(imarr[2],imarr[1]);
				trioSDoG.im[2].OP_subtract(imarr[3],imarr[2]);
				trioGradR.im[0].gradR_2(trioGradA.im[0], imarr[1]);
				trioGradR.im[1].gradR_2(trioGradA.im[1], imarr[2]);
				trioGradR.im[2].gradR_2(trioGradA.im[2], imarr[3]);
			}
			else
			{
				trioSDoG.im[1].swap(trioSDoG.im[0]);
				trioSDoG.im[2].swap(trioSDoG.im[1]);
				trioSDoG.im[2].OP_subtract(imarr[it+3],imarr[it+2]);

				trioGradR.im[1].swap(trioGradR.im[0]);
				trioGradR.im[2].swap(trioGradR.im[1]);
				trioGradA.im[1].swap(trioGradA.im[0]);
				trioGradA.im[2].swap(trioGradA.im[1]);
				trioGradR.im[2].gradR_2(trioGradA.im[2], imarr[it+3]);
			}
			calcPuntInt_act();
			calcDirecc_act();
			calcDescr_act();
		}
		if (this->imarr[nie-1].lx<=minX || this->imarr[nie-1].ly<=minY)
			break;
		nuevaOctava();
	}
	borraOctava(); // Faltaba agregar este cambio que se hizo en SDoG
	return true;
}

L_ParamManagerLocal *L_SHarris_SIFT::pideParams()
{
	return &paramSHARRISSIFT;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////
// Harris piramide
//////////////////////

L_PHarris_SIFT::L_PHarris_SIFT() : paramPHARSIFT("PHARRISSIFT",5)
{
	tipoP = PHARRIS;
	tipoD = SIFT;
	factorS=1.5;
	alfa=0.04;
	val_minP_DoG=0.01;
	val_minP_Harris=1.0e-7;
	nobj=-1;
	minX=30;
	minY=30; // de L_GenDescrPuntInt
	nPisos=0; // Pisos de la pirámide de gauss
	porcCurvMax=0.94;
	s0_def=1.0;
	s1_def=1.6;
	s1_def_chico=1.4;
	calces_maxDifIntens = 0.06;
	usa_parabol = 2;
	ultId = 0;
	// Parametros de la class_name padre
	paramPHARSIFT.addFrom("minX",&minX);
	paramPHARSIFT.addFrom("minY",&minY);
	paramPHARSIFT.addFrom("calces_maxDist",&calces_maxDist);
	paramPHARSIFT.addFrom("calces_maxRatio",&calces_maxRatio);
	paramPHARSIFT.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramPHARSIFT.addFrom("calces_probCorrecto",&calces_probCorrecto);
	paramPHARSIFT.addFrom("calces_maxDifIntens",&calces_maxDifIntens);
	// Parametros de esta class_name
	paramPHARSIFT.addFrom("factorS",&factorS);
	paramPHARSIFT.addFrom("alfa",&alfa);
	paramPHARSIFT.addFrom("val_minP_DoG",&val_minP_DoG);
	paramPHARSIFT.addFrom("val_minP_Harris",&val_minP_Harris);
	paramPHARSIFT.addFrom("porcCurvMax", &porcCurvMax);
	paramPHARSIFT.addFrom("s0_def",&s0_def);
	paramPHARSIFT.addFrom("s1_def",&s1_def);
	paramPHARSIFT.addFrom("s1_def_chico",&s1_def_chico);
	paramPHARSIFT.addFrom("usa_parabol",&usa_parabol);
	trioPDoG.f01=&genPiram.porcSub;
	trioPDoG.f12=&genPiram.porcSub;

	genPiram.paramPSS.cambiaNombreClase("PSS_PHS");
	sgen.paramSIFT.cambiaNombreClase("SIFT_PHS");
	sdirL.paramSIFTDL.cambiaNombreClase("SIFTDL_PHS");

	paramPHARSIFT.addChildren(&genPiram.paramPSS);
	paramPHARSIFT.addChildren(&sgen.paramSIFT);
	paramPHARSIFT.addChildren(&sdirL.paramSIFTDL);
}


bool L_PHarris_SIFT::calcPiramides(L_ImageGrayDouble &im, double s0)
{
	int it=0;

	genPiram.fijaImagenOriginal(im, s0);
	while (genPiram.nuevaImagen() == true)
	{
		pirLaplace[it]=genPiram.imantG;
		if (pirLaplace[it].lx<10 || pirLaplace[it].ly<10)
			break;
		pirLaplace[it]-=genPiram.imactG;
		reltam[it]=genPiram.reltam;
		pirGauss[it]=genPiram.imactG;
		pirGradR[it].gradR_2(pirGradA[it], genPiram.imantG);
		pirHarris[it].harris(genPiram.imantG,  sigm[0], factorS*sigm[0], alfa);

		it++;
		if (it>=L_PDS_NPIR)
			break;
	}
	nPisos=it;
	return true;
}

bool L_PHarris_SIFT::calcRegiones(L_ImageGrayDouble &im)
{
	int it, i, j;
	L_PuntInt pin;
	pin.imorig = &im;
	pin.imorig_lx = im.lx;
	pin.imorig_ly = im.ly;
	pin.hasDepth = 0;
	pin.tipoP = tipoP;
	L_Matrix xmax(3,1), grad, gradtr;
	L_Matrix hess, hessXY, hessinv;
	double trazaXY,detXY,val;
	double dx, dy, ds;
	xmax(0,0) = 0; xmax(1,0) = 0; xmax(2,0) = 0;

	for (it=1; it<nPisos-1; it++)
	{
		trioPDoG.im[0]=&pirLaplace[it-1];
		trioPDoG.im[1]=&pirLaplace[it];
		trioPDoG.im[2]=&pirLaplace[it+1];
		trioPHarris.im[0] = &pirHarris[it-1];
		trioPHarris.im[1] = &pirHarris[it];
		trioPHarris.im[2] = &pirHarris[it+1];

		pin.pir_n=it;
		pin.sigma0=sigm[it-1];
		pin.pir_reltam=reltam[it];

		for (j=4; j<pirLaplace[it].ly-4; j++)
		{
			for (i=4; i<pirLaplace[it].lx-4; i++)
			{
				if (pirHarris[it].pix(i,j)>=val_minP_Harris && pirHarris[it].esMax2D(i,j) && trioPDoG.esMax1D_1_abs(i,j) && L_ABS(pirLaplace[it].pix(i,j)) > val_minP_DoG && pirHarris[it].esMax2D_anillo(pirHarris[it].pix(i,j)*porcCurvMax, i, j, 5))
				{
					pin.pir_x=i;
					pin.pir_y=j;

					switch (usa_parabol)
					{
					case 0:
						trioPDoG.ajustaTaylor1((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPHarris.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break;
					case 1:
						trioPDoG.ajustaParabol((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPHarris.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break; // Error: faltaba esto
					case 2:
						dx = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i-1,j), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i+1,j));
						dy = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i,j-1), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i,j+1));
						ds = L_INTERP_MAX_PARABOL(trioPDoG.pixTrio(0,i,j), trioPDoG.pixTrio(1,i,j), trioPDoG.pixTrio(2,i,j));
						pin.pir_x += dx;
						pin.pir_y += dy;
						xmax(2,0) = ds;
						val = trioPHarris.im[1]->pix(i,j);
						break;
					}

					pin.x0=pin.pir_x/(double)reltam[it];
					pin.y0=pin.pir_y/(double)reltam[it];
					pin.intens = fabs(val);
					pin.signo = val > 0 ? 1 : -1;

					pin.sigma0=(genPiram.dupl ? genPiram.s1 : genPiram.s1_ND) *pow((double)genPiram.porcSub,(double)(-(it+ds)));

					punL.push_back(pin);
				}
			}
		}
	}
	return true;
}

bool L_PHarris_SIFT::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n;
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, pirGradR[nim], pirGradA[nim], true, sigm[0], ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_PHarris_SIFT::calcDescriptores()
{
	L_DireccNodo *p;
	int nim;
	int n=0;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n;  // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sgen.L_Descriptor::operator=( p->c );

		sgen.nobj=nobj;
		sgen.refGen=this;
		sgen.tipoP = tipoP;
		sgen.tipoD = tipoD;

		sgen.id = ultId++;
		if (sgen.calcDescriptor(pirGradR[nim], pirGradA[nim], true))
			desL.push_back(sgen);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}

bool L_PHarris_SIFT::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	int it;
	genPiram.s1=s1_def;
	genPiram.s1_ND=s1_def_chico;
	if (s0<0.1)
		s0=s0_def;

	if (genPiram.dupl)
		for (it=0; it<L_PDS_NPIR; it++)
			sigm[it]=genPiram.s1*pow((double)genPiram.porcSub,(double)(-it));
	else
		for (it=0; it<L_PDS_NPIR; it++)
			sigm[it]=genPiram.s1_ND*pow((double)genPiram.porcSub,(double)(-it));
	calcPiramides(im, s0);
	calcRegiones(im);
	calcDirecciones();
	calcDescriptores();
	return true;
}





////////////////////////////////////////////////////////////////////
// DoG piramide con descriptor SURF gaussiano

L_PDoG_SURFG::L_PDoG_SURFG() : paramPDOGSURFG("PDOGSURFG",5)
{
	tipoP = PDOG;
	tipoD = SIFT;
	r_lin=10;
	val_minP_DoG=0.01/255;
	usa_parabol=2;
	nobj=-1;
	minX=30;
	minY=30; // de L_GenDescrPuntInt
	nPisos=0; // Pisos de la pirámide de gauss
	porcCurvMax=0.94;
	s0_def=0.3;
	s1_def=1.26;
	s1_def_chico=1.26;
	calces_maxDifIntens=0.06;
	ultId = 0;
	// Parametros de la class_name padre
	paramPDOGSURFG.addFrom("minX",&minX);
	paramPDOGSURFG.addFrom("minY",&minY);
	paramPDOGSURFG.addFrom("calces_maxDist",&calces_maxDist);
	paramPDOGSURFG.addFrom("calces_maxRatio",&calces_maxRatio);
	paramPDOGSURFG.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramPDOGSURFG.addFrom("calces_probCorrecto",&calces_probCorrecto);
	// Parametros de esta class_name
	paramPDOGSURFG.addFrom("r_lin",&r_lin);
	paramPDOGSURFG.addFrom("val_minP_DoG",&val_minP_DoG);
	paramPDOGSURFG.addFrom("usa_parabol",&usa_parabol);
	paramPDOGSURFG.addFrom("porcCurvMax", &porcCurvMax);
	paramPDOGSURFG.addFrom("s0_def",&s0_def);
	paramPDOGSURFG.addFrom("s1_def",&s1_def);
	paramPDOGSURFG.addFrom("s1_def_chico",&s1_def_chico);
	trioPDoG.f01=&genPiram.porcSub;
	trioPDoG.f12=&genPiram.porcSub;

	paramPDOGSURFG.addChildren(&genPiram.paramPSS);
	paramPDOGSURFG.addChildren(&sgen.paramSURFG);
	paramPDOGSURFG.addChildren(&sdirL.paramSIFTDL);
}

bool L_PDoG_SURFG::calcPiramides(L_ImageGrayDouble &im, double s0)
{
	int it=0;

	genPiram.fijaImagenOriginal(im, s0);
	while (genPiram.nuevaImagen() == true)
	{
		pirLaplace[it]=genPiram.imantG;
		if (pirLaplace[it].lx<10 || pirLaplace[it].ly<10)
			break;
		pirLaplace[it]-=genPiram.imactG;
		reltam[it]=genPiram.reltam;
		pirGauss[it]=genPiram.imactG;
		pirGradDer[it].gradRightUp_2(pirGradArr[it], genPiram.imantG);
		pirGradR[it].gradR_2(pirGradA[it], genPiram.imantG);

		it++;
		if (it>=L_PDS_NPIR)
			break;
	}
	nPisos=it;
	return true;
}

bool L_PDoG_SURFG::calcRegiones(L_ImageGrayDouble &im)
{
	int it, i, j;
	L_PuntInt pin;
	pin.imorig = &im;
	pin.imorig_lx = im.lx;
	pin.imorig_ly = im.ly;
	pin.hasDepth = 0;
	pin.tipoP = tipoP;

	L_Matrix xmax, grad, gradtr;
	L_Matrix hess, hessXY, hessinv;
	double trazaXY,detXY,val;
	double dx, dy, ds;
	
	xmax.reallocate(3,1);
	xmax(0,0)=0; xmax(1,0)=0; xmax(2,0)=0;


	for (it=1; it<nPisos-1; it++)
	{
		trioPDoG.im[0]=&pirLaplace[it-1];
		trioPDoG.im[1]=&pirLaplace[it];
		trioPDoG.im[2]=&pirLaplace[it+1];

		pin.pir_n=it;
		pin.sigma0=sigm[it-1];
		pin.pir_reltam=reltam[it];

		for (j=4; j<pirLaplace[it].ly-4; j++)
		{
			for (i=4; i<pirLaplace[it].lx-4; i++)
			{
				if (
					(trioPDoG.im[1]->pix(i,j) > val_minP_DoG && L_PTrioIm_esMax3D(trioPDoG,i,j) && trioPDoG.im[1]->esMax2D_anillo(trioPDoG.im[1]->pix(i,j)*porcCurvMax, i, j, 3))
					||
					(trioPDoG.im[1]->pix(i,j) < -val_minP_DoG && L_PTrioIm_esMin3D(trioPDoG,i,j) && trioPDoG.im[1]->esMin2D_anillo(trioPDoG.im[1]->pix(i,j)*porcCurvMax, i, j, 3))
					)
				{
					pin.pir_x=i;
					pin.pir_y=j;
					pin.x0=i/(double)reltam[it];
					pin.y0=j/(double)reltam[it];
					pin.intens = fabs(trioPDoG.im[1]->pix(i,j));
					pin.signo = trioPDoG.im[1]->pix(i,j) > 0 ? 1 : -1;

					switch (usa_parabol)
					{
					case 0:
						trioPDoG.ajustaTaylor1((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break;
					case 1:
						trioPDoG.ajustaParabol((int)pin.pir_x,(int)pin.pir_y,grad,hess);
						detXY=hess(0,0)*hess(1,1)-hess(0,1)*hess(1,0);
						trazaXY=hess(0,0)+hess(1,1);
						gradtr.transpOf(grad);
						hessinv=hess; hessinv.invertMe();
						ds = xmax(2,0);
						val=trioPDoG.interpola(pin,xmax,grad,gradtr,hess,hessinv); // modifica pin.x0,pin.y0,pin.pir_n
						break; // Error: faltaba esto
					case 2:
						dx = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i-1,j), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i+1,j));
						dy = L_INTERP_MAX_PARABOL(trioPDoG.im[1]->pix(i,j-1), trioPDoG.im[1]->pix(i,j), trioPDoG.im[1]->pix(i,j+1));
						ds = L_INTERP_MAX_PARABOL(trioPDoG.pixTrio(0,i,j), trioPDoG.pixTrio(1,i,j), trioPDoG.pixTrio(2,i,j));
						pin.pir_x += dx;
						pin.pir_y += dy;
						xmax(2,0) = ds;
						val = trioPDoG.im[1]->pix(i,j);
						break;
					}

					pin.x0=pin.pir_x/(double)reltam[it];
					pin.y0=pin.pir_y/(double)reltam[it];
					pin.intens = fabs(val);
					pin.signo = val > 0 ? 1 : -1;

					pin.sigma0=(genPiram.dupl ? genPiram.s1 : genPiram.s1_ND) *pow((double)genPiram.porcSub,(double)(-(it+ds)));

					punL.push_back(pin);
				}
			}
		}
	}
	return true;
}

bool L_PDoG_SURFG::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n;
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, pirGradR[nim], pirGradA[nim], true, sigm[0], ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_PDoG_SURFG::calcDescriptores()
{
	L_DireccNodo *p;
	int nim;
	int n=0;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n;  // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sgen.L_Descriptor::operator=( p->c );
		sgen.id = ultId++;

		sgen.nobj=nobj;
		sgen.refGen=this;
		sgen.tipoP = tipoP;
		sgen.tipoD = tipoD;

		//if (sgen.calcDescriptor(pirGradR[nim], pirGradA[nim], true))
		//	desL.push_back(sgen);
		if (sgen.calcDescriptor(pirGradDer[nim], pirGradArr[nim], false))
			desL.push_back(sgen);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}

bool L_PDoG_SURFG::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	int i;
	genPiram.s1=s1_def;
	genPiram.s1_ND=s1_def_chico;
	if (s0<0.1)
		s0=s0_def;

	if (genPiram.dupl)
		for (i=0; i<L_PDS_NPIR; i++)
			sigm[i]=genPiram.s1*pow((double)genPiram.porcSub,(double)(-i));
	else
		for (i=0; i<L_PDS_NPIR; i++)
			sigm[i]=genPiram.s1_ND*pow((double)genPiram.porcSub,(double)(-i));
	calcPiramides(im, s0);
	calcRegiones(im);
	calcDirecciones();
	calcDescriptores();
	return true;
}

//////////////////////////////////////
// Descriptor SURF


L_SURF::L_SURF() : paramSURF("SURF",5)
{
	tipoP = P_SURF;
	tipoD = D_SURF;
	tam_0=9; deltaTam_0=3;
	porcSubm01=L_Frac(1,2);
	porcSubm12=L_Frac(1,2);
	trio.f01=&porcSubm01;
	trio.f12=&porcSubm12;
	Dcalculado=false; desSURF.tipoP=P_SURF; desSURF.tipoD=D_SURF;
	calcOri_precalc=false;
	val_minP_HoG=0.01;
	porcCurvMax=0.95;
	dupl=false;
	nrx=4;
	//nry=4;
	lx=5;
	//ly=5;
	calces_maxDifIntens = 0.1;
	ultId = 0;
	// Parametros de la class_name padre
	paramSURF.addFrom("minX",&minX);
	paramSURF.addFrom("minY",&minY);
	paramSURF.addFrom("calces_maxDist",&calces_maxDist);
	paramSURF.addFrom("calces_maxRatio",&calces_maxRatio);
	paramSURF.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramSURF.addFrom("calces_probCorrecto",&calces_probCorrecto);
	paramSURF.addFrom("calces_maxDifIntens",&calces_maxDifIntens);
	// Parametros de esta class_name
	paramSURF.addFrom("val_minP_HoG",&val_minP_HoG);
	paramSURF.addFrom("porcCurvMax",&porcCurvMax);
	paramSURF.addFrom("dupl",&dupl);
	paramSURF.addFrom("nrx",&nrx);
	paramSURF.addFrom("lx",&lx);
}

long L_SURF::genDxx(int esc, double x1, double y1, double x2, double y2)
{
	int i1, i2, j1, j2;
	//     i1   i2  0  -i2  -i1
	//  j1 *******************
	//     * P  *   N   * P  *
	//  j2 *******************
	//
	i1=L_ROUND(x1*tam[esc]); // negativo
	i2=L_ROUND(x2*tam[esc]); // negativo
	j1=L_ROUND(y1*tam[esc]); // negativo
	j2=L_ROUND(y2*tam[esc]); // positivo
	Dxx_area1_p[esc][0]=i1;
	Dxx_area1_p[esc][1]=j1;
	Dxx_area1_p[esc][2]=i2;
	Dxx_area1_p[esc][3]=j2;
	Dxx_area2_n[esc][0]=i2;
	Dxx_area2_n[esc][1]=j1;
	Dxx_area2_n[esc][2]=-i2;
	Dxx_area2_n[esc][3]=j2;
	Dxx_area3_p[esc][0]=-i2;
	Dxx_area3_p[esc][1]=j1;
	Dxx_area3_p[esc][2]=-i1;
	Dxx_area3_p[esc][3]=j2;
	return areaRectA(Dxx_area1_p[esc])+areaRectA(Dxx_area2_n[esc])+areaRectA(Dxx_area3_p[esc]);
}

long L_SURF::genDxy(int esc, double x1, double y1, double x2, double y2)
{
	int i1, i2, j1, j2;
	//     i1   i2  0  -i2  -i1
	//  j1 ******       ******
	//     * P  *       * N  *
	//  j2 ******       ******
	//
	//  0
	//
	// -j2 ******       ******
	//     * N  *       * P  *
	// -j1 ******       ******
	//
	i1=L_ROUND(x1*tam[esc]); // negativo
	i2=L_ROUND(x2*tam[esc]); // negativo
	j1=L_ROUND(y1*tam[esc]); // negativo
	j2=L_ROUND(y2*tam[esc]); // positivo
	Dxy_area1_p[esc][0]=i1;
	Dxy_area1_p[esc][1]=j1;
	Dxy_area1_p[esc][2]=i2;
	Dxy_area1_p[esc][3]=j2;
	Dxy_area2_n[esc][0]=-i2;
	Dxy_area2_n[esc][1]=j1;
	Dxy_area2_n[esc][2]=-i1;
	Dxy_area2_n[esc][3]=j2;
	Dxy_area3_n[esc][0]=i1;
	Dxy_area3_n[esc][1]=-j2;
	Dxy_area3_n[esc][2]=i2;
	Dxy_area3_n[esc][3]=-j1;
	Dxy_area4_p[esc][0]=-i2;
	Dxy_area4_p[esc][1]=-j2;
	Dxy_area4_p[esc][2]=-i1;
	Dxy_area4_p[esc][3]=-j1;
	return areaRectA(Dxy_area1_p[esc])+areaRectA(Dxy_area2_n[esc])+areaRectA(Dxy_area3_n[esc])+areaRectA(Dxy_area4_p[esc]);
}

long L_SURF::genGrX(int esc, double x1, double y1, double x2, double y2)
{
		int i1, i2, j1, j2;
	//     i1   i2  0  -i2  -i1
	//  j1 ******       ******
	//     * N  *       * P  *
	//  j2 ******       ******
	//
	i1=L_ROUND(x1*tam[esc]); // negativo
	i2=L_ROUND(x2*tam[esc]); // negativo
	j1=L_ROUND(y1*tam[esc]); // negativo
	j2=L_ROUND(y2*tam[esc]); // positivo
	GrX_area1_n[esc][0]=i1;
	GrX_area1_n[esc][1]=j1;
	GrX_area1_n[esc][2]=i2;
	GrX_area1_n[esc][3]=j2;
	GrX_area2_p[esc][0]=-i2;
	GrX_area2_p[esc][1]=j1;
	GrX_area2_p[esc][2]=-i1;
	GrX_area2_p[esc][3]=j2;
	return areaRectA(GrX_area1_n[esc])+areaRectA(GrX_area2_p[esc]);
}


void L_SURF::calcD_pn()
{
	int i;
	int deltaTam=deltaTam_0;
	tam[0]=tam_0;

	for (i=0; i<L_PDS_NPIR; i++)
	{
		if (i!=0)
		{
			tam[i]=tam[i-1]+deltaTam;
			deltaTam*=2;
			// es decir, tam[i] = tam_0 + deltaTam_0 * ( pow(2.0,i) - 1);
		}
		//printf("%d: %d\n", i, tam[i]);
		tam2[i]=tam[i]*(double)tam[i];
		tam4[i]=tam2[i]*tam2[i];
		sigma[i]=tam[i]*1.2/9.0;
		genDxx(i, -1, -5.0/9, -3.0/9, 5.0/9);
		genDxy(i, -7.0/9, -7.0/9, -1.0/9, -1.0/9);
		genGrX(i, -7.0/9, -8.0/9, 0, 8.0/9);
		//genGrX(i, -1.2*4.0, -1.2*2.0, 0, 1.2*2.0);
	}
	Dcalculado=true;
}

//#define SURF_TEST_PLZ
bool L_SURF::calcPiramides()
{
	int it, i, j;
#if defined(L_SU_B)
	#error #define conflicts
#endif
#define L_SU_B 9 //5

#ifdef SURF_TEST_PLZ
	L_ImageGrayDouble im;
	char str[100];
	original.writeBMP("orig.bmp");
#endif

	for (it=0; it<L_PDS_NPIR; it++)
	{
		pirHess[it].reallocate((int)(imIntegral.lx/sigma[it]+0.5),(int)(imIntegral.ly/sigma[it]+0.5));
		if (pirHess[it].lx<10 || pirHess[it].ly<10)
			break;
		pirGradR[it].reallocate((int)(imIntegral.lx/sigma[it]+0.5),(int)(imIntegral.ly/sigma[it]+0.5));
		pirGradA[it].reallocate((int)(imIntegral.lx/sigma[it]+0.5),(int)(imIntegral.ly/sigma[it]+0.5));
		for (j=L_SU_B; j<pirHess[it].ly-L_SU_B; j++)
		{
			for (i=L_SU_B; i<pirHess[it].lx-L_SU_B; i++)
				//pirHess[it].pix[i][j]=(double)calcDetH(it,(int)(i*sigma[it]+0.5),(int)(j*sigma[it]+0.5));
				L_SURF_calcDetH(pirHess[it].pix(i,j),imIntegral,it,(int)(i*sigma[it]+0.5),(int)(j*sigma[it]+0.5));
		}
		for (j=0; j<pirHess[it].ly; j++)
			for (i=0; i<L_SU_B; i++)
				pirHess[it].pix(i,j)=pirHess[it].pix(pirGradR[it].lx-i-1,j)=0;
		for (j=0; j<(int)(L_SU_B); j++)
			for (i=0; i<pirHess[it].lx; i++)
				pirHess[it].pix(i,j)=pirHess[it].pix(i,pirGradR[it].ly-j-1)=0;

#ifdef SURF_TEST_PLZ
		im=pirHess[it];
		//im.normalizeHistogram();
		im.multiply_to_each_element(10);
		im.sum_to_each_element(0.5);
		sprintf(str, "arch%02d.bmp", it);
		im.writeBMP(str);
#endif
	}
	nPisos=it;
#undef L_SU_B
	return true;
}

bool L_SURF::calcRegiones(L_ImageGrayDouble &im)
{
	int it, i, j;
	int bor;
	double ds;
	L_PuntInt p;
	p.imorig=&im;
	p.imorig_lx=im.lx;
	p.imorig_ly=im.ly;
	p.hasDepth = 0;
	p.tipoP=P_SURF;
	for (it=1; it<nPisos-1; it++)
	{
		p.pir_n=it;
		p.sigma0=sigma[it];
		p.pir_reltam=L_Frac(pirHess[it].lx,imIntegral.lx);
		trio.im[0]=&pirHess[it-1];
		trio.im[1]=&pirHess[it];
		trio.im[2]=&pirHess[it+1];

		// Aca es util el tener dos cambios de escala distintos en las piramides
		porcSubm01=L_Frac(tam[it-1], tam[it]);
		porcSubm12=L_Frac(tam[it], tam[it+1]);
		bor=3;

		for (j=bor; j<pirHess[it].ly-bor; j++)
		{
			for (i=bor; i<pirHess[it].lx-bor; i++)
			{
				if (
					(trio.im[1]->pix(i,j)>val_minP_HoG && L_PTrioIm_esMax3D(trio,i,j) && trio.im[1]->esMax2D_anillo(trio.im[1]->pix(i,j)*porcCurvMax,i,j,3))
					||
					(trio.im[1]->pix(i,j)<-val_minP_HoG && L_PTrioIm_esMin3D(trio,i,j) && trio.im[1]->esMin2D_anillo(trio.im[1]->pix(i,j)*porcCurvMax,i,j,3))
					)
				{
					p.pir_x = i + L_INTERP_MAX_PARABOL(trio.im[1]->pix(i-1,j), trio.im[1]->pix(i,j), trio.im[1]->pix(i+1,j));
					p.pir_y = j + L_INTERP_MAX_PARABOL(trio.im[1]->pix(i,j-1), trio.im[1]->pix(i,j), trio.im[1]->pix(i,j+1));
					ds = L_INTERP_MAX_PARABOL(trio.pixTrio(0,i,j), trio.pixTrio(1,i,j), trio.pixTrio(2,i,j));
					p.x0=p.pir_x*sigma[it];
					p.y0=p.pir_y*sigma[it];
					p.sigma0 = ( tam_0 + deltaTam_0 * ( pow(2.0,it+ds) - 1) )*1.2/9.0;
					p.intens = fabs(trio.im[1]->pix(i,j));
					p.signo = trio.im[1]->pix(i,j) > 0 ? 1 : -1;
					punL.push_back(p);
				}
			}
		}
	}
	return true;
}

bool L_SURF::vrevisaAreasIguales(int n, va_list ppp)
{
	int *a;
	long tam;
	int i;
	if (n<=0)
		return false;
	a=va_arg(ppp, int *);
	tam=(a[2]+1-a[0])*(long)(a[3]+1-a[1]);
	for (i=1; i<n; i++)
	{
		a=va_arg(ppp, int *);
		if ((a[2]+1-a[0])*(long)(a[3]+1-a[1]) != tam)
			return false;
	}
	return true;
}

double L_SURF::calcOrientacion(int esc, double i, double j, double sigma0) // Ya esta probado
{
	//  In order to be invariant to rotation, we identify a reproducible orientation for the
	//interest points. For that purpose, we first calculate the Haar-wavelet responses
	//in x and y direction, shown in Fig. 2, and this in a circular neighbourhood of
	//radius 6s around the interest point, with s the scale at which the interest point
	//was detected. Also the sampling step is scale dependent and chosen to be s. In
	//keeping with the rest, also the wavelet responses are computed at that current
	//scale s. Accordingly, at high scales the size of the wavelets is big. Therefore, we
	//use again integral images for fast filtering. Only six operations are needed to
	//compute the response in x or y direction at any scale. The side length of the
	//wavelets is 4s.
	//  Once the wavelet responses are calculated and weighted with a Gaussian (s =
	//2.5s) centered at the interest point, the responses are represented as vectors in a
	//space with the horizontal response strength along the abscissa and the vertical
	//response strength along the ordinate. The dominant orientation is estimated by
	//calculating the sum of all responses within a sliding orientation window covering
	//an angle of pi/3 . The horizontal and vertical responses within the window are
	//summed. The two summed responses then yield a new vector. The longest such
	//vector lends its orientation to the interest point. The size of the sliding window
	//is a parameter, which has been chosen experimentally. Small sizes fire on single
	//dominating wavelet responses, large sizes yield maxima in vector length that are
	//not outspoken. Both result in an unstable orientation of the interest region. Note
	//the U-SURF skips this step.

	// Version lenta pero segura

	double histOriX[12], histOriY[12], valGaus, grx, gry, a, ang;
	double sX, sY, sum, sumaMax=0, sumaMin=1e10, radio, gaU;
	int k, u, v, u0, v0, uX, vX, s0 = (int)(2.5*sigma0+0.5);

	u0=(int)(i*sigma[esc]+0.5);
	v0=(int)(j*sigma[esc]+0.5);
	radio = (6*1.42+2.5)*sigma[esc];

	if (u0 < radio || v0 < radio || u0 > imIntegral.lx-radio || v0 > imIntegral.ly-radio)
		return -1000;

	if (gaus2_5.size()==0)
		L_calcGausHNorm(gaus2_5, 2.5, 13);

	for (k=0; k<12; k++)
	{
		histOriX[k]=0;
		histOriY[k]=0;
	}


	for (u=-6; u<=6; u++)
	{
		uX = (int)(u0+u*sigma0+0.5);
		gaU = gaus2_5[u+6];
		for (v=-6; v<=6; v++)
		{
			vX = (int)(v0+v*sigma0+0.5);
			valGaus = gaU * gaus2_5[v+6];
			L_SURF_calcGrDerV_VG(grx, imIntegral, s0, uX, vX, valGaus);
			L_SURF_calcGrArrV_VG(gry, imIntegral, s0, uX, vX, valGaus);
			ang=atan2(gry, grx);
			a=(ang+2*M_PI)*12/(2*M_PI);
			k = int(a); // k < a, a > k
			histOriX[k%12]+=grx*(1-(a-k));
			histOriY[k%12]+=gry*(1-(a-k));
			histOriX[(k+1)%12]+=grx*(a-k);
			histOriY[(k+1)%12]+=gry*(a-k);
		}
	}

	ang=0;
	sX=histOriX[0]+histOriX[1]+histOriX[2];
	sY=histOriY[0]+histOriY[1]+histOriY[2];
	for (k=0; k<12; k++)
	{
		sum = sqrt(sX*sX + sY*sY);
		if (sum > sumaMax || sum < sumaMin)
		{
			if (sum > sumaMax)
			{
				sumaMax = sum;
				ang = atan2(sY, sX);
			}
			if (sum < sumaMin)
				sumaMin = sum;
		}
		sX+= histOriX[(k+3)%12] - histOriX[k%12];
		sY+= histOriY[(k+3)%12] - histOriY[k%12];
	}
	if (sumaMax < 2.0 * sumaMin)
		ang=-1000;
	return ang;
}

bool L_SURF::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	L_Direcc dir;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		dir = ppi->c;
		dir.ang=calcOrientacion(ppi->c.pir_n, ppi->c.pir_x, ppi->c.pir_y, ppi->c.sigma0);
		if (dir.ang<-100)
			continue;
		dirL.push_back(dir);
	}
	return true;

}

bool L_SURF::calcVectorDescriptor(double *vector, int esc, double i, double j, double sigma0, double ang, int &nCompVect)
{
	//  For the extraction of the descriptor, the first step consists of constructing a
	//square region centered around the interest point, and oriented along the orientation
	//selected in the previous section. For the upright version, this transformation
	//is not necessary. The size of this window is 20s. Examples of such square regions
	//are illustrated in Fig. 2.
	//  The region is split up regularly into smaller 4 × 4 square sub-regions. This
	//keeps important spatial information in. For each sub-region, we compute a few
	//simple features at 5×5 regularly spaced sample points. For reasons of simplicity,
	//we call dx the Haar wavelet response in horizontal direction and dy the Haar
	//wavelet response in vertical direction (filter size 2s). ”Horizontal” and ”vertical”
	//here is defined in relation to the selected interest point orientation. To increase
	//the robustness towards geometric deformations and localisation errors, the responses
	//dx and dy are first weighted with a Gaussian (sigma = 3.3s) centered at the
	//interest point.
	//Then, the wavelet responses dx and dy are summed up over each subregion
	//and form a first set of entries to the feature vector. In order to bring in information
	//about the polarity of the intensity changes, we also extract the sum
	//of the absolute values of the responses, |dx| and |dy|. Hence, each sub-region
	//has a four-dimensional descriptor vector v for its underlying intensity structure
	//v = (S dx, S dy, S |dx|, S |dy|). This results in a descriptor vector for all 4×4
	//sub-regions of length 64. The wavelet responses are invariant to a bias in illumination
	//(offset). Invariance to contrast (a scale factor) is achieved by turning the
	//descriptor into a unit vector.

	double r11, r12, r21, r22, valGaus, valGu, grx, gry, du, dv;
	int u, v, rx, ry, x, y, x0, y0, s0 = (int)(2*sigma0);
	double mitad, radio;

	x0=(int)(i*sigma[esc]);
	y0=(int)(j*sigma[esc]);

	radio = sigma0*(nrx*lx+1)*0.5*1.42;

	if (x0 < radio || y0 < radio || x0 >= imIntegral.lx-radio || y0 >= imIntegral.ly-radio)
		return false;

	nCompVect = 4*nrx*nrx; // 4*nrx*nry

	for (u=0; u<nCompVect; u++)
		vector[u]=0;

	// Se tiene xIm, yIm, xRot, yRot.
	// xRot = xIm (q=0) , xRot = -yIm (q=90)
	// yRot = yIm (q=0),  yRot = xIm (q=90)

	// NO, ACA ESTA EL ERROR DEL SURF
	// hay que pasar de xRot a xIm, o sea
	// xIm = xRot (q=0)  xIm = yRot  (q=90)
	// yIm = yRot (q=0)  yIm = -xRot (q=90)

	r11=cos(ang); r12=-sin(ang);
	r21=-r12;     r22=r11;


	if (gaus3_3.size() == 0)
		L_calcGausHNorm(gaus3_3,3.3, nrx*lx);

	mitad = (nrx*lx-1) / 2.0;
	for (u=0; u<nrx*lx; u++)
	{
		rx = u/lx; // Numero de region en x
		valGu = gaus3_3[u];
		for (v=0; v<nrx*lx; v++)
		{
			// Esta vuelta tiene que volar...
			ry = v/lx; // Numero de region en y
			valGaus = valGu*gaus3_3[v];

			x = (int)(sigma0*(r11*(u-mitad)+r21*(v-mitad))+0.5) + x0;
			y = (int)(sigma0*(r12*(u-mitad)+r22*(v-mitad))+0.5) + y0;

			L_SURF_calcGrDerV_VG(grx, imIntegral, s0, x, y, valGaus);
			L_SURF_calcGrArrV_VG(gry, imIntegral, s0, x, y, -valGaus);

			du = r11*grx + r12*gry;
			dv = r21*grx + r22*gry;
			vector[0 + 4*(rx+nrx*ry)]+=du;
			vector[1 + 4*(rx+nrx*ry)]+=dv;
			if (du<0)
				du=-du;
			if (dv<0)
				dv=-dv;
			vector[2 + 4*(rx+nrx*ry)]+=du;
			vector[3 + 4*(rx+nrx*ry)]+=dv;
		}
	}
	du=0;
	for (u=0; u<nCompVect; u++)
		du+=vector[u]*vector[u];
	du=1/sqrt(du);
	for (u=0; u<nCompVect; u++)
		vector[u]*=du;
	return true;
}


bool L_SURF::calcDescriptores()
{
	L_DireccNodo *p;
	int n=0;

	desSURF.nt=-1; // Sera posteriormente calculado

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		desSURF = p->c;

		desSURF.nobj=nobj;
		desSURF.refGen=this;
		desSURF.tipoP = tipoP;
		desSURF.tipoD = tipoD;

		desSURF.id = ultId++;
		if (false == calcVectorDescriptor(desSURF.vector, desSURF.pir_n, desSURF.pir_x, desSURF.pir_y, desSURF.sigma0, desSURF.ang, desSURF.nt))
			continue;
		if (dupl)
		{
			desSURF.x0/=2;
			desSURF.y0/=2;
			desSURF.sigma0/=2;
		}
		desL.push_back(desSURF);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}


bool L_SURF::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	if (dupl)
	{
		L_ImageGrayDouble im2;
		im2.doubleResolution(im);
#ifdef SURF_TEST_PLZ
		original=im2;
#endif
		imIntegral.integralOf(im2);
	}
	else
	{
#ifdef SURF_TEST_PLZ
		original=im;
#endif
		imIntegral.integralOf(im);
	}
	if (!Dcalculado)
		calcD_pn();
	calcPiramides();
	calcRegiones(im);
	calcDirecciones();
	calcDescriptores();

	return true;
}

void L_SURF::swap(L_SURF &other)
{
	other.L_GenDescrPuntInt::swap(*this);
	other.L_SURF_POD::swap(*this);
	//
	paramSURF.swap(other.paramSURF);
	imIntegral.swap(other.imIntegral);
	original.swap(other.original);
	trio.swap(other.trio);
	//
	for (int i=0; i<L_PDS_NPIR; i++)
		pirHess[i].swap(other.pirHess[i]);
	for (int i=0; i<L_PDS_NPIR; i++)
		pirGradR[i].swap(other.pirGradR[i]);
	for (int i=0; i<L_PDS_NPIR; i++)
		pirGradA[i].swap(other.pirGradA[i]);
	//
	other.desSURF.swap(desSURF);
	other.gaus2_5.swap(gaus2_5);
	other.gaus3_3.swap(gaus3_3);
}


bool L_Harris::calcRegiones(L_ImageGrayDouble &im)
{
	gradR.gradR_2(gradAng, im);
	imHarris.harris(im, 0.5, factorS*0.5, alfa);

	int i, j;
	double dx, dy;
	L_PuntInt p;
	p.hasDepth = false;
	p.tipoP = tipoP;

	for (j=1; j<imHarris.ly-1; j++)
	{
		for (i=1; i<imHarris.lx-1; i++)
		{
			if (imHarris.esMax2D(i, j, 0.95) && imHarris.pix(i,j) > val_minP_Harris)
			{
				dx = L_INTERP_MAX_PARABOL(imHarris.pix(i-1,j), imHarris.pix(i,j), imHarris.pix(i+1,j));
				dy = L_INTERP_MAX_PARABOL(imHarris.pix(i,j-1), imHarris.pix(i,j), imHarris.pix(i,j+1));
				p.x0=i+dx;
				p.y0=j+dy;
				p.sigma0 = 1;
				p.intens=fabs(imHarris.pix(i,j));
				p.signo = imHarris.pix(i,j) > 0 ? 1 : -1;
				punL.push_back(p);
			}
		}
	}
	return true;
}

bool L_Harris::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	L_DireccNodo *di;
	L_Descriptor d;


	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
		sdirL.calcDirecciones((int)ppi->c.x0, (int)ppi->c.y0, gradR, gradAng, true, 0.5, ppi);
	for (di=sdirL.lis.root; di!=NULL; di=di->sig)
	{
		d = di->c;
		d.nobj = nobj;
		d.tipoP = tipoP;
		d.tipoD = tipoD;
		d.refGen = this;
		d.id = ultId++;
		desL.push_back(d);
	}
	punL.clear();
	sdirL.lis.clear();
	return true;
}

bool L_Harris::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	calcRegiones(im);
	calcDirecciones();
	return true;
}


L_GenDescrExterno::L_GenDescrExterno(const char *lineaComando, const char *nomArchImagen, const char *nomArchDescr)
{
	tipoP = P_EXT;
	tipoD = D_EXT;
	int largoIm;
	ultId = 0;
	largoIm=(int)strlen(nomArchImagen);
	if (lineaComando != NULL)
		_lineaComando=lineaComando;
	else
		_lineaComando.resize(0);
	_nomArchImagen = nomArchImagen;
	_nomArchDescr = nomArchDescr;
	formMiko=false;
	_BMP=_PGM=false;
	if (largoIm >= 4 && _nomArchImagen[largoIm-4]=='.' &&
		(_nomArchImagen[largoIm-3] == 'B' || _nomArchImagen[largoIm-3]=='b') &&
		(_nomArchImagen[largoIm-2] == 'M' || _nomArchImagen[largoIm-2]=='m') &&
		(_nomArchImagen[largoIm-1] == 'P' || _nomArchImagen[largoIm-1]=='p')
	)
		_BMP=true;
	else if (largoIm >= 4 && _nomArchImagen[largoIm-4]=='.' &&
		(_nomArchImagen[largoIm-3] == 'P' || _nomArchImagen[largoIm-3]=='p') &&
		(_nomArchImagen[largoIm-2] == 'G' || _nomArchImagen[largoIm-2]=='g') &&
		(_nomArchImagen[largoIm-1] == 'M' || _nomArchImagen[largoIm-1]=='m')
	)
		_PGM=true;
}

L_GenDescrExterno::L_GenDescrExterno(const char *nomArchImagen)
{
	L_FileName ar, ard;

	ar.asign(nomArchImagen);
	ard.asign(nomArchImagen);
	ard.name = ard.name + "_cal";
	ard.ext = "txt";

	_lineaComando.resize(0);
	_nomArchImagen = ar.c_str();
	_nomArchDescr = ard.c_str();

	_BMP = ar.hasExtension("bmp");
	_PGM = ar.hasExtension("pgm");

	formMiko=false;
}

bool L_GenDescrExterno::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	L_ImageGrayUchar imBN;
	imBN=im;
	if (_PGM)
	{
		if (!imBN.grabaPGM(_nomArchImagen.c_str()))
		{
			printf("No se pudo crear archivo \"%s\"\n", _nomArchImagen.c_str());
			return false;
		}
	}
	else if (_BMP)
	{
		if (!imBN.writeBMP(_nomArchImagen.c_str()))
		{
			printf("No se pudo crear archivo \"%s\"\n", _nomArchImagen.c_str());
			return false;
		}
	}
	else
	{
		printf("Archivo \"%s\" con formato desconocido\n", _nomArchImagen.c_str());
		return false;
	}

	if (_lineaComando.size() != 0 && _lineaComando[0] != 0)
		if( system(_lineaComando.c_str()) ) {} // Se llama al programa externo, uso de if para eliminar warning

	if (formMiko)
		 return false;
	return desL.leeArchivoPato(_nomArchDescr.c_str(), im, this, nobj);
}

L_MinuciasExterno::L_MinuciasExterno(const char *linComando, const char *nomArchivoImagen, const char *nomArchivoMinucias):paramMinExt("MINUC_EXT",2)
{
	tipoP = MINUCIA;
	tipoD = SIFT;

	
	if (linComando == NULL && nomArchivoImagen == NULL && nomArchivoMinucias == NULL)
	{
		manejarArchivos = false; // Solo se usan los campos de orientacion y frecuencia
		lineaComando.clear();
		nomArchImagenTmp.clear();
		nomArchMinucias.clear();
		return;
	}

	manejarArchivos = true;
	lineaComando = linComando;
	nomArchImagenTmp = nomArchivoImagen;
	nomArchMinucias = nomArchivoMinucias;

	L_FileName ar(nomArchivoImagen);

	usaBMP = ar.hasExtension("bmp");
	usaPGM = ar.hasExtension("pgm");

	formMiko=false;
	agregaEspejo=true;
	usaPiramide=false;
	muestraProceso=false;
	minX = 0; // No tiene sentido aca
	minY = 0; // No tiene sentido aca

	// Parametros de la class_name padre
	paramMinExt.addFrom("minX",&minX);
	paramMinExt.addFrom("minY",&minY);
	paramMinExt.addFrom("calces_maxDist",&calces_maxDist);
	paramMinExt.addFrom("calces_maxRatio",&calces_maxRatio);
	paramMinExt.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramMinExt.addFrom("calces_probCorrecto",&calces_probCorrecto);
	// Parametros de esta class_name
	paramMinExt.addFrom("usaPiramide",&usaPiramide);
	paramMinExt.addFrom("muestraProceso",&muestraProceso);

	sgen.paramSIFT.cambiaNombreClase("SIFT_MIN");
	sdirL.paramSIFTDL.cambiaNombreClase("SIFTDL_MIN");

	paramMinExt.addChildren(&sgen.paramSIFT);
	paramMinExt.addChildren(&sdirL.paramSIFTDL);
}



L_MinuciasExterno::L_MinuciasExterno(const char *nomArchivoMinucias):paramMinExt("MINUC_EXT",2)
{
	tipoP = MINUCIA;
	tipoD = SIFT;

	throw_L_ArgException_if(nomArchivoMinucias == NULL, "L_MinuciasExterno()");

	manejarArchivos = true;
	fijaArchivoMinucias(nomArchivoMinucias);

	usaBMP = false; // No usar archivo BMP de transferencia
	usaPGM = false; // No usar archivo PGM de transferencia

	formMiko=false;
	agregaEspejo=true;
	usaPiramide=false;
	muestraProceso=false;
	minX = 0; // No tiene sentido aca
	minY = 0; // No tiene sentido aca
	// Parametros de la class_name padre
	paramMinExt.addFrom("minX",&minX);
	paramMinExt.addFrom("minY",&minY);
	paramMinExt.addFrom("calces_maxDist",&calces_maxDist);
	paramMinExt.addFrom("calces_maxRatio",&calces_maxRatio);
	paramMinExt.addFrom("calces_numCalcesPorDescr",&calces_numCalcesPorDescr);
	paramMinExt.addFrom("calces_probCorrecto",&calces_probCorrecto);
	// Parametros de esta class_name
	paramMinExt.addFrom("usaPiramide",&usaPiramide);
	paramMinExt.addFrom("muestraProceso",&muestraProceso);

	sgen.paramSIFT.cambiaNombreClase("SIFT_MIN");
	sdirL.paramSIFTDL.cambiaNombreClase("SIFTDL_MIN");

	paramMinExt.addChildren(&sgen.paramSIFT);
	paramMinExt.addChildren(&sdirL.paramSIFTDL);
}

void L_MinuciasExterno::fijaArchivoMinucias(const char *nomArchivoMinucias)
{
	L_FileName ar;

	// Solo se read el archivo de texto indicado con minucias
	ar = nomArchivoMinucias;
	ar.ext = "bmp";

	lineaComando = "";
	nomArchImagenTmp = ar.dirnameext(); // .bmp
	nomArchMinucias = nomArchivoMinucias; // .txt
}

bool L_MinuciasExterno::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	L_ImageGrayUchar imBN;
	L_ImageGrayDouble gradR, gradAng;
	imBN=im;

	throw_L_ArgException_if(manejarArchivos == false || lineaComando[0] == 0 && nomArchImagenTmp[0] == 0 && nomArchMinucias[0] == 0, "L_MinuciasExterno::calcDescrPuntInt()");

	if (lineaComando[0] != 0) // Linea de comando definida
	{
		// Grabar la imagen en el name de archivo indicado
		if (usaPGM && !imBN.grabaPGM(nomArchImagenTmp.c_str()))
		{
			printf("No se pudo crear archivo \"%s\"\n", nomArchImagenTmp.c_str());
			return false;
		}
		if (usaBMP && !imBN.writeBMP(nomArchImagenTmp.c_str()))
		{
			printf("No se pudo crear archivo \"%s\"\n", nomArchImagenTmp.c_str());
			return false;
		}
		if (muestraProceso)
			printf("  >> Llamando al programa externo...\n");
		if( system(lineaComando.c_str()) ) {} // Se llama al programa externo, uso de if para eliminar warning
	}

	if (muestraProceso)
		printf("  >> Leyendo minucias...\n");

	// Lectura de puntos de interes = minucias
	if (formMiko)
		 return false;
	if (dirL.leeArchivoPato_minu(nomArchMinucias.c_str(), im, punL, agregaEspejo) == false)
		return false;

	///// Calculo descriptores
	L_DireccNodo *dir;
	sgen.nobj=nobj;
	sgen.refGen=this;
	sgen.tipoP = tipoP;
	sgen.tipoD = tipoD;

	if (dirL.root==NULL)
		return true;
	if (usaPiramide)
	{
		L_PSS piramide;
		if (muestraProceso)
			printf("  >> Generando piramide y descriptores huella...\n");
		piramide.fijaImagenOriginal(im, s0);
		double sigmaInic=3.0;
		piramide.dupl=false;
		while(piramide.nuevaImagen() == true)
		{
			if (sigmaInic > dirL.root->c.sigma0)
			{
				gradR.gradR_2(gradAng, piramide.imactG);
				for (dir=dirL.root; dir!=NULL; dir=dir->sig)
				{
					sgen.L_Descriptor::operator=( dir->c ); // Copia tipoP, intens

					sgen.nobj=nobj;
					sgen.refGen=this;
					sgen.tipoP = tipoP;
					sgen.tipoD = tipoD;					//sgen.tipoP=FVS;

					sgen.pir_x=sgen.x0*(double)piramide.reltam;
					sgen.pir_y=sgen.y0*(double)piramide.reltam;
					sgen.ang=dir->c.ang;
					sgen.id = ultId++;
					sgen.calcDescriptor(gradR, gradAng, true);
					desL.push_back(sgen);
				}
				break;
			}
			sigmaInic/=(double)piramide.porcSub;
		}
	}
	else
	{
		if (muestraProceso)
			printf("  >> Generando descriptores huella...\n");
		gradR.gradR_2(gradAng, im);
		for (dir=dirL.root; dir!=NULL; dir=dir->sig)
		{
			sgen.L_Descriptor::operator=( dir->c );
			sgen.tipoP=tipoP;
			sgen.tipoD=DESNADA;
			sgen.pir_x=sgen.x0;
			sgen.pir_y=sgen.y0;
			sgen.ang=dir->c.ang;
			sgen.calcDescriptor(gradR, gradAng, true);
			desL.push_back(sgen);
		}
	}
	if (muestraProceso)
		printf("  >> Descriptores de minucias generados\n");
	return true;
}

bool L_MinuciasExterno::comparaCamposOrientacionFrec(double &difOri, double &difFrec, double &corrFrec, L_TransfAfinPI2D &tr, L_MinuciasExterno &other, double frecMin, double frecMax, bool usa_v, int *vx, int *vy)
{
	L_FiltroDesPolig marco(4);
	long n=0;
	double uD, vD;
	int u, v;
	double ori1_2, dori, f1, f2;
	double sumaDifAbsOri=0;
	double sumaFrec1=0;
	double sumaCuadFrec1=0;
	double sumaFrec2=0;
	double sumaCuadFrec2=0;
	double sumaFrec1Frec2=0;
	double sumaDifAbsFrec=0;
	double pf1, pf1_2, pf2, pf2_2, pf1f2;

	// No hay imagenes, resultado por defecto
	if (imFrec.data() == NULL || imOri.data() == NULL || other.imFrec.data() == NULL || other.imOri.data() == NULL)
	{
		difOri = 0; // No hay imagenes, llenar con ceros
		difFrec = 0;
		corrFrec = 0;
		return false;
	}

	int i, j;
	double lx, ly;
	if (usa_v)
	{
		marco.agregaVertice(vx[0],vy[0]);
		marco.agregaVertice(vx[1],vy[1]);
		marco.agregaVertice(vx[2],vy[2]);
		marco.agregaVertice(vx[3],vy[3]);
	}
	lx = imOri.lx;
	if (imFrec.lx < lx)
		lx=imFrec.lx;
	if (other.imOri.lx < lx)
		lx=other.imOri.lx;
	if (other.imFrec.lx < lx)
		lx=other.imFrec.lx;
	ly = imOri.ly;
	if (imFrec.ly < ly)
		ly=imFrec.ly;
	if (other.imOri.ly < ly)
		ly=other.imOri.ly;
	if (other.imFrec.ly < ly)
		ly=other.imFrec.ly;
	if (tr.m11*tr.m22 - tr.m12*tr.m21 == 0) // Transformacion absurda
		return false;
	if (tr.m21==0 && tr.m11==0) // Transformacion absurda
		return false;
	//Angulos actuales: atan2(-y, x);
	ori1_2 = -atan2(tr.m21, tr.m11); // Rotacion que apply la transformacion a las flechas de la imagen1 para dejarlas compatibles con la imagen2
	//ori1_2 = atan2(tr.m21, tr.m11); // Rotacion que apply la transformacion a las flechas de la imagen1 para dejarlas compatibles con la imagen2
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			if (imFrec.pix(i,j) < frecMin || imFrec.pix(i,j) > frecMax)
				continue;
			if (usa_v && marco.criterioEliminacion(i,j)==true)
				continue;
			uD = tr.m11*i + tr.m12*j + tr.tx;
			vD = tr.m21*i + tr.m22*j + tr.ty;
			u = (int) uD;
			v = (int) vD;
			if (u<0 || u>=lx || v<0 || v>=ly)
				continue;
			if (other.imFrec.pix(u,v) < frecMin || other.imFrec.pix(u,v) > frecMax)
				continue;
			if (usa_v && marco.criterioEliminacion(u,v)==true)
				continue;
			// El punto cumple los criterios para ser considerado en el cálculo
			dori = imOri.pix(i,j) + ori1_2 - other.imOri.pix(u,v); 
			while (dori>M_PI/2)
				dori-=M_PI;
			while (dori<-M_PI/2)
				dori+=M_PI;
			sumaDifAbsOri+=L_ABS(dori);
			f1 = imFrec.pix(i,j);
			f2 = other.imFrec.pix(u,v);
			sumaDifAbsFrec+=L_ABS(f1-f2);
			sumaFrec1+=f1;
			sumaFrec2+=f2;
			sumaCuadFrec1+=f1*f1;
			sumaCuadFrec2+=f2*f2;
			sumaFrec1Frec2+=f1*f2;
			n++;
		}
	}
	if (n==0)
		return false;
	difOri=sumaDifAbsOri/n;
	difFrec=sumaDifAbsFrec/n;
	pf1=sumaFrec1/n;
	pf1_2=sumaCuadFrec1/n;
	pf2=sumaFrec2/n;
	pf2_2=sumaCuadFrec2/n;
	pf1f2=sumaFrec1Frec2/n;
	corrFrec=(pf1f2 - pf1*pf2)/sqrt( (pf1_2 - pf1*pf1) * (pf2_2 - pf2*pf2) );
	return true;
}

void L_MinuciasExterno::grabaFlechasMinucias(const char *nomArch, L_ImageGrayDouble &im)
{
	L_ImageRGBUchar imSalida;
	L_ShapeArray lins;
	L_uchar R, G, B;
	int i, j;
	
	if (imFrec.data() == NULL || imOri.data() == NULL)
		return; // Nada que dibujar

	imSalida=im;
	for (j=0; j<imSalida.ly; j+=20)
	{
		for (i=0; i<imSalida.lx; i+=20)
		{
			if (imFrec.pix(i,j)!=0)
			{
				lins.genRandomColor(R, G, B);
				lins.drawArrow(i, j, 1/imFrec.pix(i,j), imOri.pix(i,j), R, G, B);
			}
		}
	}
	imSalida.genDrawing(lins);
	imSalida.writeBMP(nomArch);
}

bool L_MinuciasExterno::buscaCamposOrientacionFrecuencia()
{
	L_FileName archOF;
	L_ImageRGBUchar im;
	int i, j;
	double a, f;
	if (nomArchImagenTmp[0] == 0)
	{
		imFrec.destroy();
		imOri.destroy();
		return false;
	}
	archOF.asign(nomArchImagenTmp.c_str());
	archOF.name = archOF.name + "_OF";

	if (im.readImage(archOF.c_str()) != 0)
	{
		imOri.reallocate(im.lx, im.ly);
		imFrec.reallocate(im.lx, im.ly);

		for (j=0; j<im.ly; j++)
		{
			for (i=0; i<im.lx; i++)
			{
				a = im.pix(i,j,0);
				imOri.pix(i,j) = (a-127.5) * 80.53/255.0;

				f = im.pix(i,j,1);
				imFrec.pix(i,j) = (L_uchar)(f * 255.0/510.0);



			}
		}
		return true;
	}
	else
		return false; // No se encontro el archivo
}

void L_MinuciasExterno::grabaCamposOrientacionFrecuencia(FILE *fp, bool compr)
{
	L_ImageGrayUchar imf, ima;
	int i, j;
	double f, a;

	fputc(imFrec.data() != NULL, fp);
	fputc(imOri.data() != NULL, fp);

	if (imFrec.data() != NULL)
	{
		imf.reallocate(imFrec.lx, imFrec.ly);

		for (j=0; j<imf.ly; j++)
		{
			for (i=0; i<imf.lx; i++)
			{
				f = imFrec.pix(i,j);
				imf.pix(i,j) = (L_uchar)(f * 255.0/510.0);
			}
		}
		imf.grabaComprGris(fp, NULL, compr ? -1 : 0); // Busca el mejor metodo
	}

	if (imOri.data() != NULL)
	{
		ima.reallocate(imFrec.lx, imFrec.ly);

		for (j=0; j<imf.ly; j++)
		{
			for (i=0; i<imf.lx; i++)
			{
				a = imOri.pix(i,j);
				ima.pix(i,j) = (L_uchar)(127.5 + 255.0/80.53 * a);
			}
		}
		ima.grabaComprGris(fp, NULL, compr ? -1 : 0); // Busca el mejor metodo
	}
}

void L_MinuciasExterno::leeCamposOrientacionFrecuencia(FILE *fp)
{
	L_ImageGrayUchar imf, ima;
	int i, j;
	double f, a;
	bool hayf, haya;

	hayf = (fgetc(fp)==1);
	haya = (fgetc(fp)==1);

	if (hayf)
	{
		imf.leeComprGris(fp); // Busca el mejor metodo
		imFrec.reallocate(imf.lx, imf.ly);

		for (j=0; j<imf.ly; j++)
		{
			for (i=0; i<imf.lx; i++)
			{
				f = imf.pix(i,j);
				imFrec.pix(i,j) = f * 510.0 / 255.0;
			}
		}
	}

	if (haya)
	{
		ima.leeComprGris(fp); // Busca el mejor metodo
		imOri.reallocate(ima.lx, ima.ly);

		for (j=0; j<ima.ly; j++)
		{
			for (i=0; i<ima.lx; i++)
			{
				a = ima.pix(i,j);
				imOri.pix(i,j) = (a-127.5) * 80.53/255.0;
			}
		}
	}
}

bool L_GenDescrArr::agregaReferencia(L_GenDescrPuntInt *gpi)
{
	if (n>=nMax)
	{
		printf("L_GenDescrArr: nMax muy chico\n");
		return false;
	}
	arr[n]=gpi;
	activo[n]=true;
	if (calces_numCalcesPorDescr < gpi->calces_numCalcesPorDescr)
		calces_numCalcesPorDescr = gpi->calces_numCalcesPorDescr;
	paramUnion.addChildren(gpi->pideParams());
	n++;
	return true;
}


bool L_GenDescrArr::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	int i;
	bool a=true;
	for (i=0; i<n; i++)
	{
		if (activo[i])
		{
			arr[i]->nobj = nobj;
			a&=arr[i]->calcDescrPuntInt(im, s0);
		}
	}
	for (i=0; i<n; i++)
	{
		if (activo[i])
			arr[i]->desL.moveListTo(desL);
	}
	return a;
}

bool L_Harris_SIFT::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	double si, sd;
	if (dupl)
	{
		imInt.doubleResolutionBilineal(im);
		si = sigma_i_dupl;
		sd = sigma_d_dupl;
		sigma0=s0*2;
	}
	else
	{
		imInt = im;
		si = sigma_i;
		sd = sigma_d;
		sigma0=s0;
	}

	imGradR.gradR_2(imGradA, imInt);
	imHarris.harris(imInt, sd, si, alfa);

	calcRegiones(im);
	calcDirecciones();
	calcDescriptores();
	return true;
}


L_Harris_SIFT::L_Harris_SIFT():paramH1SIFT("H1SIFT",2)
{
	tipoP = HARRIS;
	tipoD = SIFT;
	nobj=-1;
	val_min_Harris = 0.0001;
	porcCurvMax = 0.95;
	alfa = 0.04;
	sigma_i = 2.0;
	sigma_d = 1.3;
	sigma_i_dupl = 4.0;
	sigma_d_dupl = 2.6;
	dupl = false;
	sMin = 4;
	ultId = 0;

	paramH1SIFT.addFrom("val_min_Harris",&val_min_Harris);
	paramH1SIFT.addFrom("porcCurvMax",&porcCurvMax);
	paramH1SIFT.addFrom("alfa",&alfa);
	paramH1SIFT.addFrom("sigma_i",&sigma_i);
	paramH1SIFT.addFrom("sigma_d",&sigma_d);
	paramH1SIFT.addFrom("sigma_i_dupl",&sigma_i_dupl);
	paramH1SIFT.addFrom("sigma_d_dupl",&sigma_d_dupl);
	paramH1SIFT.addFrom("dupl",&dupl);
	paramH1SIFT.addFrom("sMin",&sMin);

	sgen.paramSIFT.cambiaNombreClase("SIFT_H1S");
	sdirL.paramSIFTDL.cambiaNombreClase("SIFTDL_H1S");

	paramH1SIFT.addChildren(&sgen.paramSIFT);
	paramH1SIFT.addChildren(&sdirL.paramSIFTDL);
}


bool L_Harris_SIFT::calcRegiones(L_ImageGrayDouble &im)
{
	int i, j;
	L_PuntInt p;
	p.imorig=&im;
	p.imorig_lx=im.lx;
	p.imorig_ly=im.ly;
	p.hasDepth = false;
	p.tipoP=HARRIS;

	for (j=4; j<imHarris.ly-4; j++)
	{
		for (i=4; i<imHarris.lx-4; i++)
		{
			if
				(
					imHarris.esMax2D(i,j) && imHarris.pix(i,j)>=val_min_Harris
					&&
					imHarris.esMax2D_anillo(imHarris.pix(i,j)*porcCurvMax, i, j, 5)
				)
			{
				if (dupl)
				{
					p.pir_x=i/2;
					p.pir_y=j/2;
				}
				else
				{
					p.pir_x=i;
					p.pir_y=j;
				}
				p.x0=i;
				p.y0=j;
				p.intens=fabs(imHarris.pix(i,j));
				p.signo = imHarris.pix(i,j) > 0 ? 1 : -1;
				punL.push_back(p);
			}
		}
	}
	return true;
}

bool L_Harris_SIFT::calcDirecciones()
{
	L_PuntIntNodo *ppi;
	int nim;
	for (ppi=punL.root; ppi!=NULL; ppi=ppi->sig)
	{
		nim=ppi->c.pir_n;
		sdirL.calcDirecciones((int)ppi->c.pir_x, (int)ppi->c.pir_y, imGradR, imGradA, true, sigma0, ppi);
		sdirL.lis.moveListTo(dirL);
	}
	return true;
}

bool L_Harris_SIFT::calcDescriptores()
{
	L_DireccNodo *p;
	int nim;
	int n=0;

	for (p=dirL.root; p!=NULL; p=p->sig)
	{
		nim=p->c.pir_n;  // nim puede ser 0, 1 o 2 (es el numero dentro del trio)
		sgen.L_Descriptor::operator=( p->c );

		sgen.nobj=nobj;
		sgen.refGen=this;
		sgen.tipoP = tipoP;
		sgen.tipoD = tipoD;

		sgen.id = ultId++;
		if (sgen.calcDescriptor(imGradR, imGradA, true))
			desL.push_back(sgen);
		n++;
	}

	dirL.clear();
	punL.clear();
	return true;
}


bool L_MSER::calcDescrPuntInt(L_ImageGrayDouble &im, double s0)
{
	return false;
}

void L_MSER::calcRegionesMser(double valInic, double delta)
{
}






/*
void L_EvalGenDescrPuntInt::evalGPI()
{
	// Idea:
	// 1º: generar puntos de interes en imagen original => conjunto A
	// 2º: generar transformacion de semejanza e imagen semejante a la original con ruido acotado
	// 3º: generar puntos en imagen semejante => conjunto B
	// 4º: evaluar repetibilidad de puntos: ver que puntos de (AuB) corresponden a (A^B)
	// 5º: evaluar repetibilidad de descriptores de (A^B) (que tan parecidos son)
	L_ImageGrayDouble im;
	L_EvalGPIInfo inf;
	double esc, rot;

	for (esc=escMin; esc<=escMax; esc+=escDelta)
	{
		for (rot=rotMin; rot<=rotMax; rot+=rotDelta)
		{
			imOrig.cambiaTamanoFiltrando(im, pow(2.0,esc));
		}
	}
}
*/


//// Interfaces sistema antiguo




#ifdef __DIFGAUS_H_ // sistema antiguo incluido
SIFT1 * L_SHarris_SIFT::calcsift(uchar ***im3D, long lx, long ly)
{
	extern int dg_inic;
	extern SIFT1_PARAM dg;
	extern long _sift1_nobj;
	SIFT1 *ptr1;
	SIFT1 **pptr1;
	int i, j, k, u;
	L_DescriptorNodo *ptr2;
	L_SHarris_SIFT gen;
	L_ImageRGBUchar im;
	L_ImageGrayDouble imx;
	im.pix=im3D;
	im.lx=lx;
	im.ly=ly;

	if (dg_inic==0)
		inic_SIFT1_PARAM();

	gen.sgen.nx=dg.dg_nx;
	gen.sgen.ny=dg.dg_ny;
	gen.sgen.lx=dg.dg_lx;
	gen.sgen.ly=dg.dg_ly;
	gen.sgen.nd=dg.dg_nd;
	gen.val_minP=dg.dg_umbmin/4000;
	gen.s1=     dg.s0;
	gen.s1_ND=     dg.s0;

	imx=im;

#if CONVSUBM==0
	gen.sMin=100;
#elif CONVSUBM==1
	gen.sMin=0.8;
#else
#error Error en flag CONVSUBM
#endif

#if DUPLINIC==0
	gen.dupl=false;
#elif DUPLINIC==1
	gen.dupl=true;
#else
#error Error en flag DUPLINIC
#endif

	imx=im;
	gen.calcDescrPuntInt(imx, 0.5); // sigma inicial entregado

	im.pix=NULL;
	ptr2=gen.desL.root;
	ptr1=NULL;
	pptr1=&ptr1;
	while (ptr2!=NULL)
	{
		*pptr1=(SIFT1 *)malloc(sizeof(SIFT1));
		(*pptr1)->hist=callocxfloat3d(dg.dg_nd, dg.dg_nx, dg.dg_ny);
		(*pptr1)->intens=ptr2->intens; // No existe en el nuevo object
		(*pptr1)->ncap=ptr2->n0;
		(*pptr1)->esc=ptr2->reltam;
		(*pptr1)->x=(long)ptr2->x;
		(*pptr1)->y=(long)ptr2->y;
		(*pptr1)->x0=ptr2->x0;
		(*pptr1)->y0=ptr2->y0;
		(*pptr1)->ang=ptr2->ang;
		(*pptr1)->nobj=_sift1_nobj;// Se redefine despues
		u=0;
		for (i=0; i<dg.dg_nd; i++)
		{
			for (j=0; j<dg.dg_nx; j++)
			{
				for (k=0; k<dg.dg_ny; k++)
				{
					(*pptr1)->hist[i][j][k]=ptr2->vector[u++];
				}
			}
		}
		pptr1=&(*pptr1)->sig;
		ptr2=ptr2->sig;
	}
	*pptr1=NULL;
	return ptr1;
}
#endif

#ifdef __DIFGAUS_H_ // sistema antiguo incluido
SIFT1 * L_SDoG_SIFT::calcsift(uchar ***im3D, long lx, long ly)
{
	extern int dg_inic;
	extern SIFT1_PARAM dg;
	extern long _sift1_nobj;
	SIFT1 *ptr1;
	SIFT1 **pptr1;
	int i, j, k, u;
	L_DescriptorNodo *ptr2;
	L_SDoG_SIFT gen;
	L_ImageRGBUchar im;
	L_ImageGrayDouble imx;
	im.pix=im3D;
	im.lx=lx;
	im.ly=ly;

	if (dg_inic==0)
		inic_SIFT1_PARAM();

	gen.sgen.nx=dg.dg_nx;
	gen.sgen.ny=dg.dg_ny;
	gen.sgen.lx=dg.dg_lx;
	gen.sgen.ly=dg.dg_ly;
	gen.sgen.nd=dg.dg_nd;
	gen.val_min=dg.dg_umbmin;
	gen.r_lin=  dg.dg_he;
	gen.s1=     dg.s0;
	gen.s1_ND=	dg.s0;

	gen.val_minP=gen.val_min/2;
	imx=im;

#if CONVSUBM==0
	gen.sMin=100;
#elif CONVSUBM==1
	gen.sMin=0.8;
#else
#error Error en flag CONVSUBM
#endif

#if DUPLINIC==0
	gen.dupl=false;
#elif DUPLINIC==1
	gen.dupl=true;
#else
#error Error en flag DUPLINIC
#endif

	imx=im;
	gen.calcDescrPuntInt(imx, 0.5); // sigma inicial entregado

	im.pix=NULL;
	ptr2=gen.desL.root;
	ptr1=NULL;
	pptr1=&ptr1;
	while (ptr2!=NULL)
	{
		*pptr1=(SIFT1 *)malloc(sizeof(SIFT1));
		(*pptr1)->hist=callocxfloat3d(dg.dg_nd, dg.dg_nx, dg.dg_ny);
		(*pptr1)->intens=ptr2->intens; // No existe en el nuevo object
		(*pptr1)->ncap=ptr2->n0;
		(*pptr1)->esc=ptr2->reltam;
		(*pptr1)->x=(long)ptr2->x;
		(*pptr1)->y=(long)ptr2->y;
		(*pptr1)->x0=ptr2->x0;
		(*pptr1)->y0=ptr2->y0;
		(*pptr1)->ang=ptr2->ang;
		(*pptr1)->nobj=_sift1_nobj;// Se redefine despues
		u=0;
		for (i=0; i<dg.dg_nd; i++)
		{
			for (j=0; j<dg.dg_nx; j++)
			{
				for (k=0; k<dg.dg_ny; k++)
				{
					(*pptr1)->hist[i][j][k]=ptr2->vector[u++];
				}
			}
		}
		pptr1=&(*pptr1)->sig;
		ptr2=ptr2->sig;
	}
	*pptr1=NULL;
	return ptr1;
}
#endif

