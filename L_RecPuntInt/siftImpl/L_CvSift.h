#ifdef __COMPAT_CVMAT__

#ifndef __L_CVSIFT_H__
#define __L_CVSIFT_H__

#include "../L_RecPuntInt/L_GenPuntInt.h"
class L_CvSift : public L_GenDescrPuntInt
{
public:
	L_CvSift() : ultId(0)
	{
		tipoP = P_CVSIFT;
		tipoD = D_CVSIFT;

	}

	bool calcPuntInt_act(); // Calcula puntos de interés para el trio actual
	bool calcDirecc_act();  // Calcula direcciones para el trio actual
	bool calcDescr_act();   // Calcula descriptores para el trio acual
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5); // La funcion util externamente para los detectores

	long ultId;
};



#endif // __L_CVSIFT_H__

#endif // __COMPAT_CVMAT__
