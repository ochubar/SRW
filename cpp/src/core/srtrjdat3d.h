/************************************************************************//**
 * File: srtrjdat3d.h
 * Description: Electron Trajectory (and relevant characteristics) calculation in 3D Magnetic Fields (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRTRJDAT3D_H
#define __SRTRJDAT3D_H

#include "srgtrjdt.h"
#include "srmagfld.h"
#include "gmmeth.h"

//*************************************************************************

class srTTrjDat3d : public srTGenTrjDat { //Not Used? (RK moved to srTGenTrjDat)

	srTMagFld3d& m_Fld3d;
	double m_Mult2ndDer;

public:

	srTTrjDat3d(srTMagFld3d& _Fld3d) : m_Fld3d(_Fld3d) {}
	srTTrjDat3d(srTMagFld3d& _Fld3d, srTEbmDat& _EbmDat) : m_Fld3d(_Fld3d), srTGenTrjDat(&_EbmDat) {}

	//void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep); //virtual
	void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, long long ns, double sStart, double sStep); //virtual

	void funcDerivF(double s, double* arr_F, double* arr_dFds)
	{
		double xd = arr_F[1], yd = arr_F[3]; //, zd = arr_F[5];
		double zd = CGenMathMeth::radicalOnePlusSmall(-(EbmDat.GammaEm2 + xd*xd + yd*yd));
		TVector3d P(arr_F[0], arr_F[2], arr_F[4]), B;
		m_Fld3d.compB(P, B);

		arr_dFds[0] = xd;
        arr_dFds[1] = m_Mult2ndDer*(yd*B.z - zd*B.y);
		arr_dFds[2] = yd;
        arr_dFds[3] = m_Mult2ndDer*(zd*B.x - xd*B.z);
		arr_dFds[4] = zd;
        //arr_dFds[5] = m_Mult2ndDer*(xd*B.y - yd*B.x);
	}
};

//*************************************************************************

#endif
