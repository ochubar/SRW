/************************************************************************//**
 * File: srcradint.cpp
 * Description: CSR calculation
 * Project: Synchrotron Radiation Workshop
 * First release: 2006
 *
 * Copyright (C) Synchrotron SOLEIL, Gif-sur-Yvette, France
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#include "srcradint.h"
#include "sroptelm.h"
#include "srmlttsk.h"
#include "srprgind.h"
#include "gmmeth.h"

//*************************************************************************

extern srTYield srYield;

//*************************************************************************

void TAuxParamForIntegCSR::setupConstParams(srTEbmDat& eBeam)
{
	pi = 3.141592653590;
	half_pi = 0.5*pi;
	sqrt_pi = sqrt(pi);
	half_sqrt_pi = 0.5*sqrt_pi;
	two_sqrt_pi = 2*sqrt_pi;
	piE2 = pi*pi;

	double PIm10e6 = pi*1.E+06;
	half_k_d_e = PIm10e6*0.80654658; // PIm10e6dEnCon
	k_d_e = 2.*half_k_d_e;

	//OC25042017 (commented-out)
	//pxx = eBeam.GammaPartDistr('x'); pzz = eBeam.GammaPartDistr('z'); //to fix: take into account cross-moments!!!
	//pxx1 = eBeam.AlphaPartDistr('x'); pzz1 = eBeam.AlphaPartDistr('z');
	//px1x1 = eBeam.BetaPartDistr('x'); pz1z1 = eBeam.BetaPartDistr('z');
	//pxz = 0; px1z = 0; pxz1 = 0; px1z1 = 0;
	//pxg = 0; pzg = 0; px1g = 0; pz1g = 0; //to fix !
	//pxs = 0; px1s = 0; pzs = 0; pz1s = 0; //to fix !
	//psg = 0; //to fix !
	//pss = eBeam.PssPartDistr_NoCross(); //to fix !
	//pgg = eBeam.PggPartDistr_NoCross(); //to fix !

	//OC25042017
	double Mxx = eBeam.Mxx, Mxx1 = eBeam.Mxxp, Mx1x1 = eBeam.Mxpxp;
	double Myy = eBeam.Mzz, Myy1 = eBeam.Mzzp, My1y1 = eBeam.Mzpzp;
	double Mxy = eBeam.Mxz, Mx1y = eBeam.Mxpz, Mxy1 = eBeam.Mxzp, Mx1y1 = eBeam.Mxpzp;
	double Mee = eBeam.Mee; //(SigmaE/E)^2
	double Mss = eBeam.Mss;
	double Mes = eBeam.Mse;
	double Mxe = eBeam.Mxe, Mx1e = eBeam.Mxpe, Mye = eBeam.Mze, My1e = eBeam.Mzpe;
	double Mxs = eBeam.Mxs, Mx1s = eBeam.Mxps, Mys = eBeam.Mzs, My1s = eBeam.Mzps;

	double MesE2 = Mes*Mes, Mxx1E2 = Mxx1*Mxx1, Myy1E2 = Myy1*Myy1, Mxy1E2 = Mxy1*Mxy1, MxsE2 = Mxs*Mxs, MxeE2 = Mxe*Mxe, Mx1yE2 = Mx1y*Mx1y, Mx1y1E2 = Mx1y1*Mx1y1;
	double Mx1sE2 = Mx1s*Mx1s, Mx1eE2 = Mx1e*Mx1e, MysE2 = Mys*Mys, MyeE2 = Mye*Mye, My1sE2 = My1s*My1s, My1eE2 = My1e*My1e, MxyE2 = Mxy*Mxy;

	double emTotE2 = Mxx*Mx1x1*Myy*My1y1*Mee*Mss 
				    -Mx1x1*Mxx*My1y1*Myy*MesE2 - Mss*Mxx1E2*My1y1*Myy*Mee - Mss*Mx1x1*Mxx*Myy1E2*Mee - Mx1x1*My1y1*Mss*Mee*MxyE2 - Mx1x1*Myy*Mss*Mee*Mxy1E2 
				    -Mx1x1*Myy*My1y1*Mee*MxsE2 - Mx1x1*Myy*My1y1*Mss*MxeE2 - Mxx*My1y1*Mss*Mee*Mx1yE2 - Mxx*Myy*Mss*Mee*Mx1y1E2 - Mxx*Myy*My1y1*Mee*Mx1sE2 
				    -Mxx*Myy*My1y1*Mss*Mx1eE2 - Mxx*Mx1x1*My1y1*Mee*MysE2 - Mxx*Mx1x1*My1y1*Mss*MyeE2 - Mxx*Mx1x1*Myy*Mee*My1sE2 - Mxx*Mx1x1*Myy*Mss*My1eE2 
				    +Mxx*Mx1x1*Myy1E2*MesE2 + Mxx*Mx1x1*MysE2*My1eE2 + Mxx*Mx1x1*MyeE2*My1sE2 + Mxx*Myy*Mx1y1E2*MesE2 + Mxx*Myy*Mx1eE2*My1sE2 + Mxx*Myy*Mx1sE2*My1eE2 + Mxx*My1y1*Mx1yE2*MesE2 
				    +Mxx*My1y1*Mx1eE2*MysE2 + Mxx*My1y1*Mx1sE2*MyeE2 + Mxx*Mee*Mx1yE2*My1sE2 + Mxx*Mee*Mx1y1E2*MysE2 + Mxx*Mee*Mx1sE2*Myy1E2 + Mxx*Mss*Mx1yE2*My1eE2 + Mxx*Mss*Mx1y1E2*MyeE2 
				    +Mxx*Mss*Mx1eE2*Myy1E2 + Mx1x1*Myy*Mxy1E2*MesE2 + Mx1x1*Myy*MxeE2*My1sE2 + Mx1x1*Myy*MxsE2*My1eE2 + Mx1x1*My1y1*MxyE2*MesE2 + Mx1x1*My1y1*MxeE2*MysE2 + Mx1x1*My1y1*MxsE2*MyeE2 
				    +Mx1x1*Mee*MxyE2*My1sE2 + Mx1x1*Mee*Mxy1E2*MysE2 + Mx1x1*Mee*MxsE2*Myy1E2 + Mx1x1*Mss*MxyE2*My1eE2 + Mx1x1*Mss*Mxy1E2*MyeE2 + Mx1x1*Mss*MxeE2*Myy1E2 + Myy*My1y1*Mxx1E2*MesE2 
				    +Myy*My1y1*MxeE2*Mx1sE2 + Myy*My1y1*MxsE2*Mx1eE2 + Myy*Mee*Mxx1E2*My1sE2 + Myy*Mee*Mxy1E2*Mx1sE2 + Myy*Mee*MxsE2*Mx1y1E2 + Myy*Mss*Mxx1E2*My1eE2 + Myy*Mss*Mxy1E2*Mx1eE2 
				    +Myy*Mss*MxeE2*Mx1y1E2 + My1y1*Mee*Mxx1E2*MysE2 + My1y1*Mee*MxyE2*Mx1sE2 + My1y1*Mee*MxsE2*Mx1yE2 + My1y1*Mss*Mxx1E2*MyeE2 + My1y1*Mss*MxyE2*Mx1e + My1y1*Mss*MxeE2*Mx1yE2 
				    +Mss*Mee*Mxx1E2*Myy1E2 + Mss*Mee*MxyE2*Mx1y1E2 + Mss*Mee*Mxy1E2*Mx1yE2 
				    -Mxx1E2*Myy1E2*MesE2 - Mxx1E2*MyeE2*My1sE2 - Mxx1E2*MysE2*My1eE2 - MxyE2*Mx1y1E2*MesE2 - MxyE2*Mx1eE2*My1sE2 - MxyE2*Mx1sE2*My1eE2 - Mxy1E2*Mx1yE2*MesE2 - Mxy1E2*Mx1eE2*MysE2 
				    -Mxy1E2*Mx1sE2*MyeE2 - MxeE2*Mx1yE2*My1sE2 - MxeE2*Mx1y1E2*MysE2 - MxeE2*Mx1sE2*Myy1E2 - MxsE2*Mx1yE2*My1eE2 - MxsE2*Mx1y1E2*MyeE2 - MxsE2*Mx1eE2*Myy1E2;
	double emxxE2 = Mss*Mx1yE2*My1eE2 + Mx1eE2*My1y1*MysE2 + Mx1x1*My1eE2*MysE2 + Mx1eE2*My1sE2*Myy - Mss*Mx1eE2*My1y1*Myy + Mx1sE2*My1eE2*Myy - Mss*Mx1x1*My1eE2*Myy 
				   +Mss*Mx1eE2*Myy1E2 + Mss*Mx1y1E2*MyeE2 + Mx1x1*My1sE2*MyeE2 + Mx1sE2*My1y1*MyeE2 - Mss*Mx1x1*My1y1*MyeE2 + Mx1yE2*My1y1*MesE2 + Mx1y1E2*Myy*MesE2 
				   -Mx1x1*My1y1*Myy*MesE2 + Mx1x1*Myy1E2*MesE2 + Mx1yE2*My1sE2*Mee - Mss*Mx1yE2*My1y1*Mee + Mx1y1E2*MysE2*Mee - Mx1x1*My1y1*MysE2*Mee 
				   -Mss*Mx1y1E2*Myy*Mee - Mx1x1*My1sE2*Myy*Mee - Mx1sE2*My1y1*Myy*Mee + Mss*Mx1x1*My1y1*Myy*Mee + Mx1sE2*Myy1E2*Mee - Mss*Mx1x1*Myy1E2*Mee;
	double emx1x1E2 = Mss*MxyE2*My1eE2 + MxeE2*My1y1*MysE2 + Mxx*My1eE2*MysE2 + MxeE2*My1sE2*Myy - Mss*MxeE2*My1y1*Myy + MxsE2*My1eE2*Myy - Mss*Mxx*My1eE2*Myy 
					 +Mss*MxeE2*Myy1E2 + Mss*Mxy1E2*MyeE2 + Mxx*My1sE2*MyeE2 + MxsE2*My1y1*MyeE2 - Mss*Mxx*My1y1*MyeE2 + MxyE2*My1y1*MesE2 + Mxy1E2*Myy*MesE2 
					 -Mxx*My1y1*Myy*MesE2 + Mxx*Myy1E2*MesE2 + MxyE2*My1sE2*Mee - Mss*MxyE2*My1y1*Mee + Mxy1E2*MysE2*Mee - Mxx*My1y1*MysE2*Mee - Mss*Mxy1E2*Myy*Mee 
					 -Mxx*My1sE2*Myy*Mee - MxsE2*My1y1*Myy*Mee + Mss*Mxx*My1y1*Myy*Mee + MxsE2*Myy1E2*Mee - Mss*Mxx*Myy1E2*Mee;
	double emyyE2 = Mss*Mx1eE2*Mxy1E2 + Mss*Mx1y1E2*MxeE2 + Mx1eE2*Mxx*My1sE2 + Mx1x1*MxeE2*My1sE2 + Mx1eE2*MxsE2*My1y1 - Mss*Mx1eE2*Mxx*My1y1 + Mx1sE2*MxeE2*My1y1 
				   -Mss*Mx1x1*MxeE2*My1y1 + Mx1x1*MxsE2*My1eE2 + Mx1sE2*Mxx*My1eE2 - Mss*Mx1x1*Mxx*My1eE2 + Mss*Mxx1E2*My1eE2 + Mx1y1E2*Mxx*MesE2 + Mx1x1*Mxy1E2*MesE2 
				   -Mx1x1*Mxx*My1y1*MesE2 + Mxx1E2*My1y1*MesE2 + Mx1y1E2*MxsE2*Mee - Mss*Mx1y1E2*Mxx*Mee + Mx1sE2*Mxy1E2*Mee - Mss*Mx1x1*Mxy1E2*Mee - Mx1x1*Mxx*My1sE2*Mee 
				   +Mxx1E2*My1sE2*Mee - Mx1x1*MxsE2*My1y1*Mee - Mx1sE2*Mxx*My1y1*Mee + Mss*Mx1x1*Mxx*My1y1*Mee - Mss*Mxx1E2*My1y1*Mee;
	double emy1y1E2 = Mss*Mx1e*MxyE2 + Mss*Mx1yE2*MxeE2 + Mx1eE2*Mxx*MysE2 + Mx1x1*MxeE2*MysE2 + Mx1eE2*MxsE2*Myy - Mss*Mx1eE2*Mxx*Myy + Mx1sE2*MxeE2*Myy 
				     -Mss*Mx1x1*MxeE2*Myy + Mx1x1*MxsE2*MyeE2 + Mx1sE2*Mxx*MyeE2 - Mss*Mx1x1*Mxx*MyeE2 + Mss*Mxx1E2*MyeE2 + Mx1yE2*Mxx*MesE2 + Mx1x1*MxyE2*MesE2 
					 -Mx1x1*Mxx*Myy*MesE2 + Mxx1E2*Myy*MesE2 + Mx1yE2*MxsE2*Mee - Mss*Mx1yE2*Mxx*Mee + Mx1sE2*MxyE2*Mee - Mss*Mx1x1*MxyE2*Mee - Mx1x1*Mxx*MysE2*Mee 
					 +Mxx1E2*MysE2*Mee - Mx1x1*MxsE2*Myy*Mee - Mx1sE2*Mxx*Myy*Mee + Mss*Mx1x1*Mxx*Myy*Mee - Mss*Mxx1E2*Myy*Mee;
	double emeeE2 = Mss*Mx1y1E2*MxyE2 + Mss*Mx1yE2*Mxy1E2 + Mx1yE2*Mxx*My1sE2 + Mx1x1*MxyE2*My1sE2 + Mx1yE2*MxsE2*My1y1 - Mss*Mx1yE2*Mxx*My1y1 + Mx1sE2*MxyE2*My1y1 
				   -Mss*Mx1x1*MxyE2*My1y1 + Mx1y1E2*Mxx*MysE2 + Mx1x1*Mxy1E2*MysE2 - Mx1x1*Mxx*My1y1*MysE2 + Mxx1E2*My1y1*MysE2 + Mx1y1E2*MxsE2*Myy - Mss*Mx1y1E2*Mxx*Myy 
				   +Mx1sE2*Mxy1E2*Myy - Mss*Mx1x1*Mxy1E2*Myy - Mx1x1*Mxx*My1sE2*Myy + Mxx1E2*My1sE2*Myy - Mx1x1*MxsE2*My1y1*Myy - Mx1sE2*Mxx*My1y1*Myy + Mss*Mx1x1*Mxx*My1y1*Myy 
				   -Mss*Mxx1E2*My1y1*Myy + Mx1x1*MxsE2*Myy1E2 + Mx1sE2*Mxx*Myy1E2 - Mss*Mx1x1*Mxx*Myy1E2 + Mss*Mxx1E2*Myy1E2;
	double emssE2 = Mx1e*MxyE2*My1y1 + Mx1yE2*MxeE2*My1y1 + Mx1yE2*Mxx*My1eE2 + Mx1x1*MxyE2*My1eE2 + Mx1eE2*Mxy1E2*Myy + Mx1y1E2*MxeE2*Myy - Mx1eE2*Mxx*My1y1*Myy 
				   -Mx1x1*MxeE2*My1y1*Myy - Mx1x1*Mxx*My1eE2*Myy + Mxx1E2*My1eE2*Myy + Mx1eE2*Mxx*Myy1E2 + Mx1x1*MxeE2*Myy1E2 + Mx1y1E2*Mxx*MyeE2 + Mx1x1*Mxy1E2*MyeE2 
				   -Mx1x1*Mxx*My1y1*MyeE2 + Mxx1E2*My1y1*MyeE2 + Mx1y1E2*MxyE2*Mee + Mx1yE2*Mxy1E2*Mee - Mx1yE2*Mxx*My1y1*Mee - Mx1x1*MxyE2*My1y1*Mee - Mx1y1E2*Mxx*Myy*Mee 
				   -Mx1x1*Mxy1E2*Myy*Mee + Mx1x1*Mxx*My1y1*Myy*Mee - Mxx1E2*My1y1*Myy*Mee - Mx1x1*Mxx*Myy1E2*Mee + Mxx1E2*Myy1E2*Mee;
	double emxx1E2 = My1eE2*MysE2 - Mss*My1eE2*Myy + My1sE2*MyeE2 - Mss*My1y1*MyeE2 - My1y1*Myy*MesE2 + Myy1E2*MesE2 - My1y1*MysE2*Mee - My1sE2*Myy*Mee + Mss*My1y1*Myy*Mee - Mss*Myy1E2*Mee;
	double emxyE2 = Mx1eE2*My1sE2 - Mss*Mx1eE2*My1y1 + Mx1sE2*My1eE2 - Mss*Mx1x1*My1eE2 + Mx1y1E2*MesE2 - Mx1x1*My1y1*MesE2 - Mss*Mx1y1E2*Mee - Mx1x1*My1sE2*Mee - Mx1sE2*My1y1*Mee + Mss*Mx1x1*My1y1*Mee;
	double emxy1E2 = Mx1eE2*MysE2 - Mss*Mx1eE2*Myy + Mx1sE2*MyeE2 - Mss*Mx1x1*MyeE2 + Mx1yE2*MesE2 - Mx1x1*Myy*MesE2 - Mss*Mx1yE2*Mee - Mx1x1*MysE2*Mee - Mx1sE2*Myy*Mee + Mss*Mx1x1*Myy*Mee;
	double emxeE2 = Mx1yE2*My1sE2 - Mss*Mx1yE2*My1y1 + Mx1y1E2*MysE2 - Mx1x1*My1y1*MysE2 - Mss*Mx1y1E2*Myy - Mx1x1*My1sE2*Myy - Mx1sE2*My1y1*Myy + Mss*Mx1x1*My1y1*Myy + Mx1sE2*Myy1E2 - Mss*Mx1x1*Myy1E2;
	double emxsE2 = Mx1yE2*My1eE2 - Mx1eE2*My1y1*Myy - Mx1x1*My1eE2*Myy + Mx1eE2*Myy1E2 + Mx1y1E2*MyeE2 - Mx1x1*My1y1*MyeE2 - Mx1yE2*My1y1*Mee - Mx1y1E2*Myy*Mee + Mx1x1*My1y1*Myy*Mee - Mx1x1*Myy1E2*Mee;
	double emx1yE2 = MxeE2*My1sE2 - Mss*MxeE2*My1y1 + MxsE2*My1eE2 - Mss*Mxx*My1eE2 + Mxy1E2*MesE2 - Mxx*My1y1*MesE2 - Mss*Mxy1E2*Mee - Mxx*My1sE2*Mee - MxsE2*My1y1*Mee + Mss*Mxx*My1y1*Mee;
	double emx1y1E2 = MxeE2*MysE2 - Mss*MxeE2*Myy + MxsE2*MyeE2 - Mss*Mxx*MyeE2 + MxyE2*MesE2 - Mxx*Myy*MesE2 - Mss*MxyE2*Mee - Mxx*MysE2*Mee - MxsE2*Myy*Mee + Mss*Mxx*Myy*Mee;
	double emx1eE2 = MxyE2*My1sE2 - Mss*MxyE2*My1y1 + Mxy1E2*MysE2 - Mxx*My1y1*MysE2 - Mss*Mxy1E2*Myy - Mxx*My1sE2*Myy - MxsE2*My1y1*Myy + Mss*Mxx*My1y1*Myy + MxsE2*Myy1E2 - Mss*Mxx*Myy1E2;
	double emx1sE2 = MxyE2*My1eE2 - MxeE2*My1y1*Myy - Mxx*My1eE2*Myy + MxeE2*Myy1E2 + Mxy1E2*MyeE2 - Mxx*My1y1*MyeE2 - MxyE2*My1y1*Mee - Mxy1E2*Myy*Mee + Mxx*My1y1*Myy*Mee - Mxx*Myy1E2*Mee;
	double emyy1E2 = Mx1eE2*MxsE2 - Mss*Mx1eE2*Mxx + Mx1sE2*MxeE2 - Mss*Mx1x1*MxeE2 - Mx1x1*Mxx*MesE2 + Mxx1E2*MesE2 - Mx1x1*MxsE2*Mee - Mx1sE2*Mxx*Mee + Mss*Mx1x1*Mxx*Mee - Mss*Mxx1E2*Mee;
	double emyeE2 = Mx1y1E2*MxsE2 - Mss*Mx1y1E2*Mxx + Mx1sE2*Mxy1E2 - Mss*Mx1x1*Mxy1E2 - Mx1x1*Mxx*My1sE2 + Mxx1E2*My1sE2 - Mx1x1*MxsE2*My1y1 - Mx1sE2*Mxx*My1y1 + Mss*Mx1x1*Mxx*My1y1 - Mss*Mxx1E2*My1y1;
	double emysE2 = Mx1eE2*Mxy1E2 + Mx1y1E2*MxeE2 - Mx1eE2*Mxx*My1y1 - Mx1x1*MxeE2*My1y1 - Mx1x1*Mxx*My1eE2 + Mxx1E2*My1eE2 - Mx1y1E2*Mxx*Mee - Mx1x1*Mxy1E2*Mee + Mx1x1*Mxx*My1y1*Mee - Mxx1E2*My1y1*Mee;
	double emy1eE2 = Mx1yE2*MxsE2 - Mss*Mx1yE2*Mxx + Mx1sE2*MxyE2 - Mss*Mx1x1*MxyE2 - Mx1x1*Mxx*MysE2 + Mxx1E2*MysE2 - Mx1x1*MxsE2*Myy - Mx1sE2*Mxx*Myy + Mss*Mx1x1*Mxx*Myy - Mss*Mxx1E2*Myy;
	double emy1sE2 = Mx1e*MxyE2 + Mx1yE2*MxeE2 - Mx1eE2*Mxx*Myy - Mx1x1*MxeE2*Myy - Mx1x1*Mxx*MyeE2 + Mxx1E2*MyeE2 - Mx1yE2*Mxx*Mee - Mx1x1*MxyE2*Mee + Mx1x1*Mxx*Myy*Mee - Mxx1E2*Myy*Mee;
	double emesE2 = Mx1y1E2*MxyE2 + Mx1yE2*Mxy1E2 - Mx1yE2*Mxx*My1y1 - Mx1x1*MxyE2*My1y1 - Mx1y1E2*Mxx*Myy - Mx1x1*Mxy1E2*Myy + Mx1x1*Mxx*My1y1*Myy - Mxx1E2*My1y1*Myy - Mx1x1*Mxx*Myy1E2 + Mxx1E2*Myy1E2;

	double inv2_emTotE2 = 0.5/emTotE2;
	pxx = emxxE2*inv2_emTotE2;
	px1x1 = emx1x1E2*inv2_emTotE2;
	pzz = emyyE2*inv2_emTotE2;
	pz1z1 = emy1y1E2*inv2_emTotE2;
	pgg = emeeE2*inv2_emTotE2;
	pss = emssE2*inv2_emTotE2;
	pxx1 = (-Mxx1*emxx1E2)*inv2_emTotE2;
	pxz = (-Mxy*emxyE2)*inv2_emTotE2;
	pxz1 = (-Mxy1*emxy1E2)*inv2_emTotE2;
	pxg = (-Mxe*emxeE2)*inv2_emTotE2;
	pxs = (-Mxs*emxsE2)*inv2_emTotE2;
	px1z = (-Mx1y*emx1yE2)*inv2_emTotE2;
	px1z1 = (-Mx1y1*emx1y1E2)*inv2_emTotE2;
	px1g = (-Mx1e*emx1eE2)*inv2_emTotE2;
	px1s = (-Mx1s*emx1sE2)*inv2_emTotE2;
	pzz1 = (-Myy1*emyy1E2)*inv2_emTotE2;
	pzg = (-Mye*emyeE2)*inv2_emTotE2;
	pzs = (-Mys*emysE2)*inv2_emTotE2;
	pz1g = (-My1e*emy1eE2)*inv2_emTotE2;
	pz1s = (-My1s*emy1sE2)*inv2_emTotE2;
	psg = (-Mes*emesE2)*inv2_emTotE2;

	constW1 = 2*(pxx*eBeam.x0 + pxx1*eBeam.dxds0 + pxz*eBeam.z0 + pxz1*eBeam.dzds0 + pxs*eBeam.sc);
	constW2 = 2*(pxx1*eBeam.x0 + px1x1*eBeam.dxds0 + px1z*eBeam.z0 + px1z1*eBeam.dzds0 + px1s*eBeam.sc);
	constW3a = 2*(pxz*eBeam.x0 + px1z*eBeam.dxds0 + pzz*eBeam.z0 + pzz1*eBeam.dzds0 + pzs*eBeam.sc);
	constW3b = 2*(pxz1*eBeam.x0 + px1z1*eBeam.dxds0 + pzz1*eBeam.z0 + pz1z1*eBeam.dzds0 + pz1s*eBeam.sc);
	constW4b = 2*(pxs*eBeam.x0 + px1s*eBeam.dxds0 + pzs*eBeam.z0 + pz1s*eBeam.dzds0 + pss*eBeam.sc);
	constW5 = 2*(pxg*eBeam.x0 + px1g*eBeam.dxds0 + pzg*eBeam.z0 + pz1g*eBeam.dzds0 + psg*eBeam.sc);
	constQ6 = -pxx*eBeam.x0*eBeam.x0 - 2.*pxx1*eBeam.x0*eBeam.dxds0 - px1x1*eBeam.dxds0*eBeam.dxds0 -2.*pxz*eBeam.x0*eBeam.z0 
		-2.*px1z*eBeam.dxds0*eBeam.z0 - pzz*eBeam.z0*eBeam.z0 - 2.*pxz1*eBeam.x0*eBeam.dzds0 - 2.*px1z1*eBeam.dxds0*eBeam.dzds0 
		-2.*pzz1*eBeam.z0*eBeam.dzds0 - pz1z1*eBeam.dzds0*eBeam.dzds0 -2.*pxs*eBeam.x0*eBeam.sc - 2.*px1s*eBeam.dxds0*eBeam.sc 
		-2.*pzs*eBeam.z0*eBeam.sc - 2.*pz1s*eBeam.dzds0*eBeam.sc - pss*eBeam.sc*eBeam.sc;

	//OC26042017 (uncommented)
	double f02 = -pxx; //OC030110
	double f02e2 = f02*f02; //OC030110
	double f12 = -((f02*px1x1 + pxx1*pxx1)/f02); //OC030110
	double k0 = -(pxx1*pxg) - f02*px1g; //OC030110
	double k1 = f02*px1z + pxx1*pxz; //OC030110
	double k3 = f02*px1z1 + pxx1*pxz1; //OC030110
	double k2 = (k0*k1)/(f02e2*f12) - (pxz*pxg)/f02 - pzg; //OC030110
	double k4 = -(k1*k3)/(f02e2*f12) - (pxz*pxz1)/f02 - pzz1; //OC030110
	double k5 = (pxz*pxg)/f02 + (k1*(pxx1*pxg + f02*px1g))/(f02e2*f12) + pzg; //OC030110
	double k6 = (pxz*pxs)/f02 - (k1*(-(f02*px1s) - pxx1*pxs))/(f02e2*f12) + pzs; //OC030110
	double f22 = -(k1*k1/f12 + f02*(pxz*pxz + f02*pzz))/f02e2; //OC030110
    double k7 = (k0*k3)/(f02e2*f12) + (k4*k5)/f22 - (pxz1*pxg)/f02 - pz1g; //OC030110
	double k8 = (k4*k6)/f22 - (pxz1*pxs)/f02 + (k3*(-(f02*px1s) - pxx1*pxs))/(f02e2*f12) - pz1s; //OC030110
	double f32 = -(k3*k3/(f02e2*f12)) - pxz1*pxz1/f02 - (k4*k4 + f22*pz1z1)/f22; //OC030110
	double k9 = (2*k2*k6)/f22 - (2*k7*k8)/f32 - (2*k0*(-(f02*px1s) - pxx1*pxs))/(f02e2*f12) - (2*pxs*pxg)/f02 - 2*psg; //OC030110
	double buf42_1 = f02*px1s + pxx1*pxs; //OC030110
	double buf42_2 = f02e2*f12*k4*k4 + f22*(k3*k3 + f02*f12*(pxz1*pxz1 + f02*pz1z1)); //OC030110
	double f42 = -(k6*k6/f22) - pxs*pxs/f02 - buf42_1*buf42_1/(f02e2*f12) + (f02e2*f12*f22*k8*k8)/buf42_2 - pss;
	double f52 = -(k0*k0/(f02e2*f12)) - k2*k2/f22 - k7*k7/f32 - k9*k9/(4.*f42) - pxg*pxg/f02 - pgg; //OC030110
	cf = sqrt(f02*f12*f22*f32*f42*f52)/(pi*pi*pi);

    //double cf_no_cross = sqrt(pxx*px1x1*pzz*pz1z1*pss*pgg)/(pi*pi*pi); //test
	//OC26042017 (commented-out)
    //cf = sqrt(pxx*px1x1*pzz*pz1z1*pss*pgg)/(pi*pi*pi); //test

	//double sig_z1 = sqrt(1./(2.*pz1z1));
	//double sig_z = sqrt(1./(2.*pzz));

	const double alpha = 1./137.0360411; // Fine-structure constant
	const double e_coulomb = 1.602189246E-19; // Charge of electron in C
	const double convPhEn = 0.80654658;
	const double auxCn = alpha*convPhEn*convPhEn*(1.E+03)/e_coulomb;
	cn = cf*sqrt(auxCn*eBeam.Current*(::fabs(eBeam.Neb - 1.))); // To be multiplied by photon energy in [eV]
}

//*************************************************************************

void srTCSR::checkInputConsistency()
{
	//throw err_no if there are problems
}

//*************************************************************************

void srTCSR::estimateAbsoluteTolerance()
{//To check this !!! 

	//to implement
/**
	double GammaE2 = (TrjDatPtr->EbmDat.Gamma)*(TrjDatPtr->EbmDat.Gamma);
	double AvgObsDist = (DistrInfoDat.yStart + DistrInfoDat.yEnd)*0.5;
	double AvgObsWaveLength = (DistrInfoDat.LambStart + DistrInfoDat.LambEnd)*0.5;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*AvgObsWaveLength : PIm10e6*1000./AvgObsWaveLength;

	double c10e9_d_Lamb = PIm10e9_d_Lamb/PI;
	double InvDiffractAngleE2 = AvgObsDist*c10e9_d_Lamb;
	double Min_GammaE2_InvDiffractAngleE2 = (GammaE2 < InvDiffractAngleE2)? GammaE2 : InvDiffractAngleE2;
	double InvAvgObsDist = 1./AvgObsDist;
	double ExpectedIntesityValues = (3.47E+12)*(TrjDatPtr->EbmDat.Current)*Min_GammaE2_InvDiffractAngleE2*InvAvgObsDist*InvAvgObsDist;

	EstimatedAbsoluteTolerance = (1.E-05)*ExpectedIntesityValues*sIntegRelPrec;
**/
}

//*************************************************************************

void srTCSR::analyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;
	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;

	//TrjDatPtr->AnalizeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
	// CW 3,4 fails to make this call. So the function is here:
	if(!m_TrjDat.HorFieldIsNotZero) FieldIsSymOverZ = 1;
	if(!m_TrjDat.VerFieldIsNotZero) FieldIsSymOverX = 1;

	if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;

	char ObsIsSymOverX = 0, ObsIsSymOverZ = 0;
	if(FieldIsSymOverX && (m_Wfr.nx > 1))
	{
		double xTol = m_Wfr.xStep*0.01; // To steer
		char TrjAngIsSmall = (::fabs(m_TrjDat.EbmDat.dxds0*(m_Wfr.yStart - m_TrjDat.EbmDat.s0)) < xTol);
        ObsIsSymOverX = TrjAngIsSmall && (::fabs(m_Wfr.GetWfrMiddleHor() - m_TrjDat.EbmDat.x0) < xTol);
        //ObsIsSymOverX = TrjAngIsSmall && (::fabs(0.5*(m_Wfr.xStart + DistrInfoDat.xEnd) - m_TrjDat.EbmDat.x0) < xTol);
	}
	if(FieldIsSymOverZ && (m_Wfr.nz > 1))
	{
		//double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
		double zTol = m_Wfr.zStep*0.01; // To steer
		char TrjAngIsSmall = (::fabs(m_TrjDat.EbmDat.dzds0*(m_Wfr.yStart - m_TrjDat.EbmDat.s0)) < zTol);
		ObsIsSymOverZ = TrjAngIsSmall && (::fabs(m_Wfr.GetWfrMiddleVer() - m_TrjDat.EbmDat.z0) < zTol);
		//ObsIsSymOverZ = TrjAngIsSmall && (::fabs(0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd) - TrjDatPtr->EbmDat.z0) < zTol);
	}
	if((!ObsIsSymOverX) && (!ObsIsSymOverZ)) return;

	char ElecBeamIsSymOverX = 0, ElecBeamIsSymOverZ = 0;
	m_TrjDat.EbmDat.analyseDistribSymmetry(ElecBeamIsSymOverX, ElecBeamIsSymOverZ);

	FinalResAreSymOverX = (FieldIsSymOverX && ObsIsSymOverX && ElecBeamIsSymOverX);
	FinalResAreSymOverZ = (FieldIsSymOverZ && ObsIsSymOverZ && ElecBeamIsSymOverZ);
}

//*************************************************************************

void srTCSR::computeTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont)
{
    srTFieldBasedArrayKeys Keys;

	if(pMagLensCont != 0)
	{
		computeOffAxisTrajArrays(FldArr, pMagLensCont);
		Keys.ZeroAllKeys();
		Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = Keys.IntBtxE2_ = Keys.IntBtzE2_ = 1;
		Keys.IntX01toS0_ = Keys.IntX02toS0_ = Keys.IntX11toS0_ = Keys.IntX12toS0_ = Keys.IntX22toS0_ = 1;
		Keys.IntZ01toS0_ = Keys.IntZ02toS0_ = Keys.IntZ11toS0_ = Keys.IntZ12toS0_ = Keys.IntZ22toS0_ = 1;
		Keys.IntX01_ = Keys.IntX02_ = Keys.IntX11_ = Keys.IntX12_ = Keys.IntX22_ = Keys.IntZ01_ = Keys.IntZ02_ = Keys.IntZ11_ = Keys.IntZ12_ = Keys.IntZ22_ = 1;
		m_TrjDat.CompTotalTrjData(Keys, FldArr);
	}
	else
	{
		Keys.ZeroAllKeys();
		Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = Keys.IntBtxE2_ = Keys.IntBtzE2_ = 1;
		Keys.X1p_ = Keys.Z1p_ = Keys.X2p_ = Keys.Z2p_ = Keys.X1_ = Keys.Z1_ = Keys.X2_ = Keys.Z2_ = 1;
		Keys.IntX01toS0_ = Keys.IntX02toS0_ = Keys.IntX11toS0_ = Keys.IntX12toS0_ = Keys.IntX22toS0_ = 1;
		Keys.IntZ01toS0_ = Keys.IntZ02toS0_ = Keys.IntZ11toS0_ = Keys.IntZ12toS0_ = Keys.IntZ22toS0_ = 1;
		Keys.IntX01_ = Keys.IntX02_ = Keys.IntX11_ = Keys.IntX12_ = Keys.IntX22_ = Keys.IntZ01_ = Keys.IntZ02_ = Keys.IntZ11_ = Keys.IntZ12_ = Keys.IntZ22_ = 1;
		m_TrjDat.CompTotalTrjData(Keys, FldArr);
	}
}

//*************************************************************************

void srTCSR::computeOffAxisTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont)
{
	if(pMagLensCont == 0) return;

    double *pX1p = FldArr.X1pArr, *pZ1p = FldArr.Z1pArr;
    double *pX2p = FldArr.X2pArr, *pZ2p = FldArr.Z2pArr;
    double *pX1 = FldArr.X1Arr, *pZ1 = FldArr.Z1Arr;
	double *pX2 = FldArr.X2Arr, *pZ2 = FldArr.Z2Arr;

	double s = FldArr.sStart;
	double sStp = FldArr.sStep;
	for(int i=0; i<FldArr.Ns; i++)
	{
        TMatrix2d Mx(1,0,1,0), Mz(1,0,1,0);
		pMagLensCont->ComputeParticlePropagMatrix(s, Mx, Mz);

		*(pX1++) = Mx.Str0.x; *(pX2++) = Mx.Str0.y;
		*(pX1p++) = Mx.Str1.x; *(pX2p++) = Mx.Str1.y;
		*(pZ1++) = Mz.Str0.x; *(pZ2++) = Mz.Str0.y;
		*(pZ1p++) = Mz.Str1.x; *(pZ2p++) = Mz.Str1.y;
		s += sStp;
	}
}

//*************************************************************************

void srTCSR::copySymEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, char SymWithRespectToXax, char SymWithRespectToZax)
{
	char ChangeSignEx = SymWithRespectToZax;
	char ChangeSignEz = SymWithRespectToXax;

	float *tOrigEx = pOrigDataEx, *tSymEx = pSymDataEx;
	float *tOrigEz = pOrigDataEz, *tSymEz = pSymDataEz;
	for(int ie=0; ie<m_Wfr.ne; ie++)
	{
		*tSymEx = *(tOrigEx++); *(tSymEx + 1) = *(tOrigEx++);
		if(ChangeSignEx) { *tSymEx = -(*tSymEx); *(tSymEx + 1) *= -1;}
		tSymEx += 2;

		*tSymEz = *(tOrigEz++); *(tSymEz + 1) = *(tOrigEz++);
		if(ChangeSignEz) { *tSymEz = -(*tSymEz); *(tSymEz + 1) *= -1;}
		tSymEz += 2;
	}
}

//*************************************************************************

void srTCSR::fillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ)
{
	//long PerX = m_Wfr.ne << 1;
	//long PerZ = PerX*m_Wfr.nx;
	long long PerX = m_Wfr.ne << 1;
	long long PerZ = PerX*m_Wfr.nx;

	char SymWithRespectToXax, SymWithRespectToZax;
	int HalfNz = m_Wfr.nz >> 1, Nz_mi_1 = m_Wfr.nz - 1;
	//int izStart = ((HalfNz << 1) == m_Wfr.nz)? HalfNz : (HalfNz + 1); //OC030110
	int HalfNx = m_Wfr.nx >> 1, Nx_mi_1 = m_Wfr.nx - 1;
	//int ixStart = ((HalfNx << 1) == m_Wfr.nx)? HalfNx : (HalfNx + 1);
	int iz, ix;

	if(FinalResAreSymOverZ)
	{
		if(FinalResAreSymOverX)
		{
			SymWithRespectToXax = 0; SymWithRespectToZax = 1;
			for(iz=0; iz<HalfNz; iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				for(ix=0; ix<HalfNx; ix++)
				{
					//long Offset = izPerZ + ix*PerX;
					long long Offset = izPerZ + ix*PerX;
					float* pOrigDataEx = m_Wfr.pBaseRadX + Offset;
					float* pOrigDataEz = m_Wfr.pBaseRadZ + Offset;
					//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					float* pSymDataEx = m_Wfr.pBaseRadX + OffsetSym;
					float* pSymDataEz = m_Wfr.pBaseRadZ + OffsetSym;
					copySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
		SymWithRespectToXax = 1; SymWithRespectToZax = 0;
		for(iz=0; iz<HalfNz; iz++)
		{
			//long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
			long long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
			for(ix=0; ix<m_Wfr.nx; ix++)
			{
				//long ixPerX = ix*PerX;
				//long Offset = izPerZ + ixPerX;
				long long ixPerX = ix*PerX;
				long long Offset = izPerZ + ixPerX;
				float* pOrigDataEx = m_Wfr.pBaseRadX + Offset;
				float* pOrigDataEz = m_Wfr.pBaseRadZ + Offset;
				//long OffsetSym = BufZ + ixPerX;
				long long OffsetSym = BufZ + ixPerX;
				float* pSymDataEx = m_Wfr.pBaseRadX + OffsetSym;
				float* pSymDataEz = m_Wfr.pBaseRadZ + OffsetSym;
				copySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
			}
		}
	}
	else if(FinalResAreSymOverX)
	{
		SymWithRespectToXax = 0; SymWithRespectToZax = 1;
		for(iz=0; iz<m_Wfr.nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			for(ix=0; ix<HalfNx; ix++)
			{
				//long Offset = izPerZ + ix*PerX;
				long long Offset = izPerZ + ix*PerX;
				float* pOrigDataEx = m_Wfr.pBaseRadX + Offset;
				float* pOrigDataEz = m_Wfr.pBaseRadZ + Offset;
				//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				float* pSymDataEx = m_Wfr.pBaseRadX + OffsetSym;
				float* pSymDataEz = m_Wfr.pBaseRadZ + OffsetSym;
				copySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
			}
		}
	}
}

//*************************************************************************

void srTCSR::setupInitialTrajArrays(srTMagFldCont* pMagLensCont)
{
    //if((pTrUnifTrjDat == 0) && (pMagLensCont == 0)) throw INCORRECT_PARAMS_SR_COMP;
    //if(pPrcPar == 0) throw INCORRECT_PARAMS_SR_COMP;
	//if(pPrcPar->RelPrecOrStep <= 0) throw INCORRECT_STEP_OR_RELATIVE_PRECISION;

	if(m_PrecParams.m_MethNo == 0)
	{// constant grid, fixed step vs s
		//m_FldArr.sStep = m_PrecParams.m_PrecPar;
		m_FldArr.sStep = m_PrecParams.m_sStep;

		m_FldArr.sStart = m_TrjDat.sStart;
		if(m_FldArr.sStart < m_PrecParams.m_sIntegStart) m_FldArr.sStart = m_PrecParams.m_sIntegStart;
		double sEnd = m_TrjDat.sStart + (m_TrjDat.LenFieldData	- 1)*(m_TrjDat.sStep);
		if(sEnd > m_PrecParams.m_sIntegEnd) sEnd = m_PrecParams.m_sIntegEnd;

		m_FldArr.Ns = int((sEnd - m_FldArr.sStart)/(m_FldArr.sStep));
        m_FldArr.Ns = ((m_FldArr.Ns >> 1) << 1) + 1; // to ensure odd

		if(m_FldArr.Ns < 5) throw TOO_LARGE_LONGITUD_INTEG_STEP;
		m_FldArr.sStep = (sEnd - m_FldArr.sStart)/(m_FldArr.Ns - 1);

		srTFieldBasedArrayKeys K;
		K.Btx_ = K.Btz_ = K.X_ = K.Z_ = K.IntBtxE2_ = K.IntBtzE2_ = 1;
        K.X1p_ = K.Z1p_ = K.X2p_ = K.Z2p_ = K.X1_ = K.Z1_ = K.X2_ = K.Z2_ = 1;
        K.IntX01_ = K.IntX02_ = K.IntX11_ = K.IntX12_ = K.IntX22_ = K.IntZ01_ = K.IntZ02_ = K.IntZ11_ = K.IntZ12_ = K.IntZ22_ = 1;
        m_FldArr.AllocateArrays(m_FldArr.Ns, K);
        computeTrajArrays(m_FldArr, pMagLensCont);

	//	AllocateCoefArraysForInteg2D(gFldArr.Ns);
	//	AllocateFuncArraysForExternInteg(gFldArr.Ns);
	}
	else if(m_PrecParams.m_MethNo == 1)
	{
        //program new method here
	}
}

//*************************************************************************

void srTCSR::performMethodDependentSetupActions()
{
	srTMagFldCont* pMagLensCont = 0;
    setupInitialTrajArrays(pMagLensCont);

	m_AuxIntPar.setupConstParams(m_TrjDat.EbmDat);
	m_AuxIntPar.allocateArrays(m_PrecParams.m_MethNo, m_FldArr.Ns);

	if(m_PrecParams.m_MC_NumMacroPart > 0) 
	{// Monte-Carlo, for tests
		m_pRadInt = new srTRadInt();
		m_pWfrSmpAux = new srTWfrSmp(m_Wfr.yStart, m_Wfr.xStart, m_Wfr.xStart + m_Wfr.xStep*(m_Wfr.nx - 1), 1, m_Wfr.zStart, m_Wfr.zStart + m_Wfr.zStep*(m_Wfr.nz - 1), 1, 0, m_Wfr.eStart, m_Wfr.eStart + m_Wfr.eStep*(m_Wfr.ne - 1), 1, "EV");
		m_pWfrAux = new srTSRWRadStructAccessData(&(m_TrjDat.EbmDat), &m_TrjDat, m_pWfrSmpAux, 0);
		
		double prec_or_step = (m_PrecParams.m_MethNo > 0)? m_PrecParams.m_PrecPar : m_PrecParams.m_sStep;
		m_pParPrecElecFldSingle = new srTParPrecElecFld(m_PrecParams.m_MethNo, prec_or_step, m_PrecParams.m_sIntegStart, m_PrecParams.m_sIntegEnd, m_PrecParams.m_NxNzOversampFact);

		m_TrjDat.DeallocateMemoryForCfs();

		m_pTrjDatAux = new srTTrjDat(m_TrjDat);
		m_pTrjDatAux->m_doNotDeleteData = true;

		srTEbmDat& e_beam = m_TrjDat.EbmDat;
		m_xcArr[0] = e_beam.x0;
		m_xcArr[1] = e_beam.dxds0;
		m_xcArr[2] = e_beam.z0;
		m_xcArr[3] = e_beam.dzds0;
		m_xcArr[4] = e_beam.Energy;
		m_xcArr[5] = e_beam.sc;
        m_sigArr[0] = e_beam.GetSigma_x();
        m_sigArr[1] = e_beam.GetSigma_dxds();
        m_sigArr[2] = e_beam.GetSigma_z();
        m_sigArr[3] = e_beam.GetSigma_dzds();
        m_sigArr[4] = e_beam.GetSigmaE_GeV();
        m_sigArr[5] = e_beam.GetSigma_s();
	}

	if(m_PrecParams.m_MethNo == 0) 
	{
		//radIntegrationManual(exzy, dEwdsAtEdges, Ew);
	}
    else if(m_PrecParams.m_MethNo == 1) 
	{
		//radIntegrationAutoUnd(exzy, dEwdsAtEdges, Ew);
		estimateAbsoluteTolerance();
	}
    else if(m_PrecParams.m_MethNo == 2) 
	{
		//radIntegrationAutoWig(exzy, dEwdsAtEdges, Ew);
		estimateAbsoluteTolerance();
	}
}

//*************************************************************************

void srTCSR::performMethodDependentFinishActions()
{
	m_AuxIntPar.disposeArrays();

	if(m_PrecParams.m_MC_NumMacroPart > 0) 
	{// Monte-Carlo, for tests
		if(m_pRadInt != 0) { delete m_pRadInt; m_pRadInt = 0;}
		if(m_pWfrSmpAux != 0) { delete m_pWfrSmpAux; m_pWfrSmpAux = 0;}
		if(m_pWfrAux != 0) { delete m_pWfrAux; m_pWfrAux = 0;}
		if(m_pParPrecElecFldSingle != 0) { delete m_pParPrecElecFldSingle; m_pParPrecElecFldSingle = 0;}
		if(m_pTrjDatAux != 0) { delete m_pTrjDatAux; m_pTrjDatAux = 0;}
	}

	if(m_PrecParams.m_MethNo == 0) 
	{
		//radIntegrationManual(exzy, dEwdsAtEdges, Ew);
	}
    else if(m_PrecParams.m_MethNo == 1) 
	{
		//radIntegrationAutoUnd(exzy, dEwdsAtEdges, Ew);
	}
    else if(m_PrecParams.m_MethNo == 2) 
	{
		//radIntegrationAutoWig(exzy, dEwdsAtEdges, Ew);
	}
}

//*************************************************************************

void srTCSR::computeElectricFieldFreqDomain()
{
	checkInputConsistency();
	performMethodDependentSetupActions();

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	analyzeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ);

	double StepE = (m_Wfr.ne > 1)? m_Wfr.eStep : 0.;
	double StepX = (m_Wfr.nx > 1)? m_Wfr.xStep : 0.;
	double StepZ = (m_Wfr.nz > 1)? m_Wfr.zStep : 0.;

	double xc = m_TrjDat.EbmDat.x0;
	double zc = m_TrjDat.EbmDat.z0;
	double xTol = StepX*0.001, zTol = StepZ*0.001; // To steer

	//long PerX = m_Wfr.ne << 1;
	//long PerZ = m_Wfr.nx*PerX;
	long long PerX = m_Wfr.ne << 1;
	long long PerZ = m_Wfr.nx*PerX;
	float *pEx0 = m_Wfr.pBaseRadX;
	float *pEz0 = m_Wfr.pBaseRadZ;

	srTEXZY instEXZY;
	instEXZY.y = m_Wfr.yStart;
	instEXZY.z = m_Wfr.zStart;
	srTEFourier Ew;

	//long TotalAmOfOutPoints = m_Wfr.nz*m_Wfr.nx*m_Wfr.ne, PointCount = 0;
	long long TotalAmOfOutPoints = m_Wfr.nz*m_Wfr.nx*m_Wfr.ne, PointCount = 0;
	if(FinalResAreSymOverX) TotalAmOfOutPoints >>= 1;
	if(FinalResAreSymOverZ) TotalAmOfOutPoints >>= 1;
	double UpdateTimeInt_s = 0.5;
    srTCompProgressIndicator CompProgressInd(TotalAmOfOutPoints, UpdateTimeInt_s);

	int result = 0;
	for(int iz=0; iz<m_Wfr.nz; iz++)
	{
        if(FinalResAreSymOverZ) { if((instEXZY.z - zc) > zTol) break;}

		//long izPerZ = iz*PerZ;
		long long izPerZ = iz*PerZ;
		instEXZY.x = m_Wfr.xStart;
		for(int ix=0; ix<m_Wfr.nx; ix++)
		{
			if(FinalResAreSymOverX) { if((instEXZY.x - xc) > xTol) break;}

			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;
			instEXZY.e = m_Wfr.eStart;
			for(int ie=0; ie<m_Wfr.ne; ie++)
			{
				genRadIntegration(instEXZY, Ew);

				//long Offset = izPerZ + ixPerX + (ie << 1);
				long long Offset = izPerZ + ixPerX + (ie << 1);
				float *pEx = pEx0 + Offset, *pEz = pEz0 + Offset;

				*pEx = (float)Ew.EwX_Re; //float(RadIntegValues->real());
				*(pEx+1) = (float)Ew.EwX_Im; //float(RadIntegValues->imag());
				*pEz = (float)Ew.EwZ_Re; //float(RadIntegValues[1].real());
				*(pEz+1) = (float)Ew.EwZ_Im; //float(RadIntegValues[1].imag());

				if(result = CompProgressInd.UpdateIndicator(PointCount++)) throw result;
				if(result = srYield.Check()) throw result;

				instEXZY.e += StepE;
			}
			instEXZY.x += StepX;
		}
		instEXZY.z += StepZ;
	}
	if(FinalResAreSymOverZ || FinalResAreSymOverX) fillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ);

	performMethodDependentFinishActions();

	srTGenOptElem GenOptElem;
	int res = 0;
	if(res = GenOptElem.ComputeRadMoments(&m_Wfr)) throw res;
}

//*************************************************************************

void srTCSR::radIntegrationResiduals(srTEXZY& exzy,	srTEFourier& Ew, srTEFourier* arr_dEwds)
{//assumes that m_FldArr.sStep, etc. parameters are already defined
	srTEFourier	EwDummy;
	const int numberOfTerms = 2; //not more than 2!
	const int npDer = 3; //can't be less than 3!
	complex<double> arrAmpX[npDer], arrAmpZ[npDer], arrArg[npDer];
	complex<double> arr_dAXds[3], arr_dAZds[3], arr_dBds[4];
	complex<double> arrAuxDer[2];
	complex<double> expB, dEwXds_L, dEwXds_R, dEwZds_L, dEwZds_R;

//Left Residual
	for(int i=0; i<npDer; i++) 
	{
        computeFuncToIntegAtOnePointOnTrj(i, exzy, EwDummy, arrAmpX[i], arrAmpZ[i], arrArg[i]);
	}
	arr_dAXds[0] = arrAmpX[0];
	arr_dAXds[1] = CGenMathMeth::differentiate(arrAmpX, m_FldArr.sStep, 0, npDer);
	arrAuxDer[0] = arr_dAXds[1]; arrAuxDer[1] = CGenMathMeth::differentiate(arrAmpX, m_FldArr.sStep, 1, npDer);
	arr_dAXds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 0, 2); //d2AXds2

	arr_dAZds[0] = arrAmpZ[0];
	arr_dAZds[1] = CGenMathMeth::differentiate(arrAmpZ, m_FldArr.sStep, 0, npDer);
	arrAuxDer[0] = arr_dAZds[1]; arrAuxDer[1] = CGenMathMeth::differentiate(arrAmpZ, m_FldArr.sStep, 1, npDer);
	arr_dAZds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 0, 2); //d2AZds2

	arr_dBds[0] = arrArg[0];
	arr_dBds[1] = CGenMathMeth::differentiate(arrArg, m_FldArr.sStep, 0, npDer);
    arrAuxDer[0] = arr_dBds[1]; arrAuxDer[1] = CGenMathMeth::differentiate(arrArg, m_FldArr.sStep, 1, npDer);
	arr_dBds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 0, 2); //d2Bds2

	complex<double> EwXL = computeResidualRow(arr_dAXds, arr_dBds, numberOfTerms);
	complex<double> EwZL = computeResidualRow(arr_dAZds, arr_dBds, numberOfTerms);

	if(arr_dEwds != 0)
	{
        expB = exp(arr_dBds[0]);
        dEwXds_L = (arr_dAXds[1] + (arr_dAXds[0]*arr_dBds[1]))*expB;
        dEwZds_L = (arr_dAZds[1] + (arr_dAZds[0]*arr_dBds[1]))*expB;
	}

//Right Residual
	int iCount = 0;
	//for(int i=(m_FldArr.Ns - npDer); i<m_FldArr.Ns; i++) 
	for(long long i=(m_FldArr.Ns - npDer); i<m_FldArr.Ns; i++) 
	{
        computeFuncToIntegAtOnePointOnTrj(i, exzy, EwDummy, arrAmpX[iCount], arrAmpZ[iCount], arrArg[iCount]);
		iCount++;
	}
	int npDer_mi_1 = npDer - 1, npDer_mi_2 = npDer - 2; //, npDer_mi_3 = npDer - 3; //OC030110

    arr_dAXds[0] = arrAmpX[npDer_mi_1];
    arr_dAXds[1] = CGenMathMeth::differentiate(arrAmpX, m_FldArr.sStep, npDer_mi_1, npDer); //dAXds
	arrAuxDer[1] = arr_dAXds[1]; arrAuxDer[0] = CGenMathMeth::differentiate(arrAmpX, m_FldArr.sStep, npDer_mi_2, npDer);
    arr_dAXds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 1, 2); //d2AXds2

    arr_dAZds[0] = arrAmpZ[npDer_mi_1];
    arr_dAZds[1] = CGenMathMeth::differentiate(arrAmpZ, m_FldArr.sStep, npDer_mi_1, npDer); //dAZds
    arrAuxDer[1] = arr_dAZds[1]; arrAuxDer[0] = CGenMathMeth::differentiate(arrAmpZ, m_FldArr.sStep, npDer_mi_2, npDer);
	arr_dAZds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 1, 2); //d2AZds2

    arr_dBds[0] = arrArg[npDer_mi_1];
    arr_dBds[1] = CGenMathMeth::differentiate(arrArg, m_FldArr.sStep, npDer_mi_1, npDer); //dBds
    arrAuxDer[1] = arr_dBds[1]; arrAuxDer[0] = CGenMathMeth::differentiate(arrArg, m_FldArr.sStep, npDer_mi_2, npDer);
	arr_dBds[2] = CGenMathMeth::differentiate(arrAuxDer, m_FldArr.sStep, 1, 2); //d2Bds2

    complex<double> EwXR = computeResidualRow(arr_dAXds, arr_dBds, numberOfTerms);
    complex<double> EwZR = computeResidualRow(arr_dAZds, arr_dBds, numberOfTerms);
	
	//complex<double> dEwX = EwXR - EwXL, dEwZ = EwZR - EwZL;
	complex<double> dEwX = EwXL - EwXR, dEwZ = EwZL - EwZR;

	Ew.EwX_Re = dEwX.real(); Ew.EwX_Im = dEwX.imag();
	Ew.EwZ_Re = dEwZ.real(); Ew.EwZ_Im = dEwZ.imag();

	if(arr_dEwds != 0)
	{
        expB = exp(arr_dBds[0]);
        dEwXds_R = (arr_dAXds[1] + (arr_dAXds[0]*arr_dBds[1]))*expB;
        dEwZds_R = (arr_dAZds[1] + (arr_dAZds[0]*arr_dBds[1]))*expB;

		srTEFourier &dEwds_L = arr_dEwds[0], &dEwds_R = arr_dEwds[1];
		dEwds_L.EwX_Re = dEwXds_L.real(); dEwds_L.EwX_Im = dEwXds_L.imag();
		dEwds_L.EwZ_Re = dEwZds_L.real(); dEwds_L.EwZ_Im = dEwZds_L.imag();
		dEwds_R.EwX_Re = dEwXds_L.real(); dEwds_R.EwX_Im = dEwXds_R.imag();
		dEwds_R.EwZ_Re = dEwZds_L.real(); dEwds_R.EwZ_Im = dEwZds_R.imag();
	}
}

//*************************************************************************

void srTCSR::radIntegrationManual(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew)
{
	//double s = m_FldArr.sStart;
	srTEFourier initDer, finDer, auxDer;
	complex<double> auxAmpX, auxAmpZ, auxArg;

	//long Ns_mi_1 = m_FldArr.Ns - 1;
	long long Ns_mi_1 = m_FldArr.Ns - 1;
	//for(long i=0; i<m_FldArr.Ns; i++)
	for(long long i=0; i<m_FldArr.Ns; i++)
	{
		char calcDeriv = 0;
		if((i == 0) || (i == Ns_mi_1)) calcDeriv = 1;

		computeFuncToIntegAtOnePointOnTrj(i, exzy, m_AuxIntPar.arrEw[i], auxAmpX, auxAmpZ, auxArg);

		//s += m_FldArr.sStep;

		//if(i == 0) initDer = calcDeriv;
		//else if(i == Ns_mi_1) finDer = calcDeriv;
	}
    integrateSimpleEwArr(m_AuxIntPar.arrEw, m_FldArr.Ns, m_FldArr.sStep, arr_dEwds, Ew);
}

//*************************************************************************

//void srTCSR::integrateSimpleEwArr(srTEFourier* arrEw, long np, double h, srTEFourier* pDer, srTEFourier& resEw)
void srTCSR::integrateSimpleEwArr(srTEFourier* arrEw, long long np, double h, srTEFourier* pDer, srTEFourier& resEw)
{//Np is assumed odd and >= 5!!!
	const double we = 7./15.;
    const double w1 = 16./15.;
    const double w2 = 14./15.;
    const double wd = 1./15.; 

	srTEFourier Sum1 = arrEw[1];
	srTEFourier Sum2(0,0,0,0);

	srTEFourier *t = arrEw + 2;
	//long nLoops = (np - 3) >> 1;
	long long nLoops = (np - 3) >> 1;
    //for(long i=0; i<nLoops; i++)
    for(long long i=0; i<nLoops; i++)
	{
		Sum2 += *(t++);
		Sum1 += *(t++);
	}

	srTEFourier *pFe1 = arrEw;
	srTEFourier *pFe2 = arrEw + (np - 1);

	Sum1 *= w1;
	Sum2 *= w2;
	Sum1 += Sum2;

	double he2_wd = h*h*wd;
	resEw = h*((we*(*pFe1 + *pFe2)) + Sum1);
	if(pDer != 0)
	{
		resEw += (he2_wd*(*pDer)) + (he2_wd*(*(pDer + 1)));
	}
}

//*************************************************************************

void srTCSR::radIntegrationAutoUnd(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew)
{

}

//*************************************************************************

void srTCSR::radIntegrationAutoWig(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew)
{

}

//*************************************************************************

//void srTCSR::computeFuncToIntegAtOnePointOnTrj(long i, srTEXZY exzy, srTEFourier& Ew, complex<double>& ampX, complex<double>& ampZ, complex<double>& arg)
void srTCSR::computeFuncToIntegAtOnePointOnTrj(long long i, srTEXZY exzy, srTEFourier& Ew, complex<double>& ampX, complex<double>& ampZ, complex<double>& arg)
{
	double s = m_FldArr.sStart + i*m_FldArr.sStep;
	double invR = 1./(exzy.y - s);

	//double *pBx = m_FldArr.BxArr + i, *pBz = m_FldArr.BzArr + i;
    double *pX = m_FldArr.XArr + i, *pZ = m_FldArr.ZArr + i;
	double *pBtx = m_FldArr.BtxArr + i, *pBtz = m_FldArr.BtzArr + i;
	double *pX1 = m_FldArr.X1Arr + i, *pZ1 = m_FldArr.Z1Arr + i, *pX1p = m_FldArr.X1pArr + i, *pZ1p = m_FldArr.Z1pArr + i;
	double *pX2 = m_FldArr.X2Arr + i, *pZ2 = m_FldArr.Z2Arr + i, *pX2p = m_FldArr.X2pArr + i, *pZ2p = m_FldArr.Z2pArr + i;
	double *pIx10 = m_FldArr.IntX01Arr + i, *pIz10 = m_FldArr.IntZ01Arr + i;
	double *pIx20 = m_FldArr.IntX02Arr + i, *pIz20 = m_FldArr.IntZ02Arr + i;
	double *pIx11 = m_FldArr.IntX11Arr + i, *pIz11 = m_FldArr.IntZ11Arr + i;
	double *pIx12 = m_FldArr.IntX12Arr + i, *pIz12 = m_FldArr.IntZ12Arr + i;
	double *pIx22 = m_FldArr.IntX22Arr + i, *pIz22 = m_FldArr.IntZ22Arr + i;

    double k = exzy.e*m_AuxIntPar.k_d_e;
	double k_d_2 = 0.5*k;

    double x1dr = (*pX1)*invR, z1dr = (*pZ1)*invR, x2dr = (*pX2)*invR, z2dr = (*pZ2)*invR;
	double i00 = (*(m_FldArr.IntBtxE2Arr + i)) + (*(m_FldArr.IntBtzE2Arr + i));
	double x_mi_xeq0 = exzy.x - *pX, z_mi_zeq0 = exzy.z - *pZ;
	double s_gamEm2 = s*m_TrjDat.EbmDat.GammaEm2;

	double phi_0 = k_d_2*(i00 + (x_mi_xeq0*x_mi_xeq0 + z_mi_zeq0*z_mi_zeq0)*invR + s_gamEm2);
	double phi_g = k*(-i00 + (x_mi_xeq0*(*pX) + z_mi_zeq0*(*pZ))*invR - s_gamEm2);
	
	//double phi_s = -k*m_TrjDat.EbmDat.sc;
	double phi_s = -k;

	double phi_x = k*((*pIx10) - x_mi_xeq0*x1dr), phi_z = k*((*pIz10) - z_mi_zeq0*z1dr);
	double phi_x1 = k*((*pIx20) - x_mi_xeq0*x2dr), phi_z1 = k*((*pIz20) - z_mi_zeq0*z2dr);
	double phi_xx = k_d_2*((*pIx11) + x1dr*(*pX1)), phi_zz = k_d_2*((*pIz11) + z1dr*(*pZ1));
	double phi_xx1 = k*((*pIx12) + x1dr*(*pX2)), phi_zz1 = k*((*pIz12) + z1dr*(*pZ2));
	double phi_x1x1 = k_d_2*((*pIx22) + x2dr*(*pX2)), phi_z1z1 = k_d_2*((*pIz22) + z2dr*(*pZ2));
	double phi_xg = -k*((*pIx10) + x1dr*(*pX)), phi_zg = -k*((*pIz10) + z1dr*(*pZ));
	double phi_x1g = -k*((*pIx20) + x2dr*(*pX)), phi_z1g = -k*((*pIz20) + z2dr*(*pZ));
	double phi_gg = k_d_2*(i00 + ((*pX)*(*pX) + (*pZ)*(*pZ))*invR);

	double invRE2 = invR*invR;
	double x_mi_xeq0_d_RE2 = x_mi_xeq0*invRE2, z_mi_zeq0_d_RE2 = z_mi_zeq0*invRE2;
	double xpeq0_d_R = (*pBtx)*invR, zpeq0_d_R = (*pBtz)*invR;

	complex<double> bufC(1., invR/k);
	complex<double> h0 = ((-x_mi_xeq0_d_RE2)*bufC) + xpeq0_d_R;
	complex<double> hg = ((-(*pX)*invRE2)*bufC) - xpeq0_d_R;
	complex<double> hx = ((x1dr*invR)*bufC) + ((*pX1p)*invR);
	complex<double> hx1 = ((x2dr*invR)*bufC) + ((*pX2p)*invR);
	complex<double> v0 = ((-z_mi_zeq0_d_RE2)*bufC) + zpeq0_d_R;
	complex<double> vg = ((-(*pZ)*invRE2)*bufC) - zpeq0_d_R;
	complex<double> vz = ((z1dr*invR)*bufC) + ((*pZ1p)*invR);
	complex<double> vz1 = ((z2dr*invR)*bufC) + ((*pZ2p)*invR);

	complex<double> sxx(-m_AuxIntPar.pxx, phi_xx), szz(-m_AuxIntPar.pzz, phi_zz);
	complex<double> sxx1(-2*m_AuxIntPar.pxx1, phi_xx1), szz1(-2*m_AuxIntPar.pzz1, phi_zz1); //OC26042017 (uncommented)
		//complex<double> sxx1(0, phi_xx1), szz1(0, phi_zz1);

	complex<double> sx1x1(-m_AuxIntPar.px1x1, phi_x1x1), sz1z1(-m_AuxIntPar.pz1z1, phi_z1z1);
	complex<double> sxg(-2*m_AuxIntPar.pxg, phi_xg), szg(-2*m_AuxIntPar.pzg, phi_zg); //OC26042017 (uncommented)
		//complex<double> sxg(0, phi_xg), szg(0, phi_zg);

	complex<double> s0 = sxx;
	complex<double> sxx1_d_s0 = sxx1/s0;
	complex<double> sx1z = (-2*m_AuxIntPar.px1z) + (m_AuxIntPar.pxz*sxx1_d_s0); //OC26042017 (uncommented)
		//complex<double> sx1z = 0;

	complex<double> sx1z1 = (-2*m_AuxIntPar.px1z1) + (m_AuxIntPar.pxz1*sxx1_d_s0); //OC26042017 (uncommented)
		//complex<double> sx1z1 = 0;

	complex<double> sx1s = (-2*m_AuxIntPar.px1s) + (m_AuxIntPar.pxs*sxx1_d_s0); //OC26042017 (uncommented)
		//complex<double> sx1s = 0;

	complex<double> s1 = sx1x1 - (0.25*sxx1*sxx1_d_s0);

	complex<double> pxz_d_s0 = m_AuxIntPar.pxz/s0, sx1z_d_s1 = sx1z/s1; //OC26042017 (uncommented)
		//complex<double> pxz_d_s0 = 0, sx1z_d_s1 = 0;

	complex<double> pxz1_d_s0 = m_AuxIntPar.pxz1/s0, sx1z1_d_s1 = sx1z1/s1; //OC26042017 (uncommented)
		//complex<double> pxz1_d_s0 = 0, sx1z1_d_s1 = 0;

	complex<double> szz1a = szz1 - ((2*m_AuxIntPar.pxz1)*pxz_d_s0) - (0.5*sx1z1*sx1z_d_s1); //OC26042017 (uncommented)
		//complex<double> szz1a = szz1;

	complex<double> szs = (-2*m_AuxIntPar.pzs) - ((2*m_AuxIntPar.pxs)*pxz_d_s0) - (0.5*sx1s*sx1z_d_s1); //OC26042017 (uncommented)
		//complex<double> szs = 0;

	complex<double> s2 = szz - (m_AuxIntPar.pxz*pxz_d_s0) - (0.25*sx1z*sx1z_d_s1); //OC26042017 (uncommented)
		//complex<double> s2 = szz;

	complex<double> szz1a_d_s2 = szz1a/s2;
	
	complex<double> sz1s = (-2*m_AuxIntPar.pz1s) - (2*m_AuxIntPar.pxs*pxz1_d_s0) - (0.5*sx1s*sx1z1_d_s1) - (0.5*szs*szz1a_d_s2); //OC26042017 (uncommented)
		//complex<double> sz1s = 0;

    complex<double> sx1g(-2*m_AuxIntPar.px1g, phi_x1g); //OC26042017 (uncommented)
		//complex<double> sx1g(0, phi_x1g);

	sx1g -= 0.5*sxg*sxx1_d_s0;

	complex<double> szga = (-0.5*sx1g*sx1z_d_s1) + (sxg*pxz_d_s0) + szg; //OC26042017 (uncommented)
		//complex<double> szga = szg;

    complex<double> sz1g(-2*m_AuxIntPar.pz1g, phi_z1g); //OC26042017 (uncommented)
		//complex<double> sz1g(0, phi_z1g);

	sz1g += (sxg*pxz1_d_s0) - (0.5*sx1g*sx1z1_d_s1) - (0.5*szga*szz1a_d_s2); //OC26042017 (uncommented)
		//sz1g += -(0.5*szga*szz1a_d_s2);

	complex<double> s3 = sz1z1 - (m_AuxIntPar.pxz1*pxz1_d_s0) - (0.25*sx1z1*sx1z1_d_s1) - (0.25*szz1a*szz1a_d_s2); //OC26042017 (uncommented)
		//complex<double> s3 = sz1z1 - (0.25*szz1a*szz1a_d_s2);

	complex<double> ssg = (-2*m_AuxIntPar.psg) + (m_AuxIntPar.pxs*sxg/s0) - (0.5*sx1s*sx1g/s1) - (0.5*szs*szga/s2) - (0.5*sz1s*sz1g/s3); //OC26042017 (uncommented)
		//complex<double> ssg = 0;
	
	complex<double> s4 = (-m_AuxIntPar.pss) - ((m_AuxIntPar.pxs*m_AuxIntPar.pxs)/s0) - (0.25*((sx1s*sx1s/s1) + (szs*szs/s2) + (sz1s*sz1s/s3))); //OC26042017 (uncommented)
		//complex<double> s4 = (-m_AuxIntPar.pss);

	complex<double> s5(-m_AuxIntPar.pgg, phi_gg);
	
	s5 -= 0.25*((sxg*sxg/s0) + (sx1g*sx1g/s1) + (szga*szga/s2) + (sz1g*sz1g/s3) + (ssg*ssg/s4)); //OC26042017 (uncommented)
		//s5 -= 0.25*((sxg*sxg/s0) + (sx1g*sx1g/s1) + (szga*szga/s2) + (sz1g*sz1g/s3));

	complex<double> a0 = hx;
	complex<double> sqrt_mi_s0 = sqrtC(-s0), sqrt_mi_s1 = sqrtC(-s1), sqrt_mi_s2 = sqrtC(-s2), sqrt_mi_s3 = sqrtC(-s3), sqrt_mi_s4 = sqrtC(-s4), sqrt_mi_s5 = sqrtC(-s5);
	complex<double> mi_s0_e32 = -s0*sqrt_mi_s0;
	complex<double> inv_mi_s0_e32 = 1./mi_s0_e32, inv_mi_s1_e32 = 1./(-s1*sqrt_mi_s1), inv_mi_s2_e32 = 1./(-s2*sqrt_mi_s2), inv_mi_s3_e32 = 1./(-s3*sqrt_mi_s3), inv_mi_s4_e32 = 1./(-s4*sqrt_mi_s4), inv_mi_s5_e32 = 1./(-s5*sqrt_mi_s5);

	complex<double> a1 = (m_AuxIntPar.half_sqrt_pi*inv_mi_s0_e32)*((a0*sxx1) - (2.*hx1*s0));
	complex<double> a0_s1_inv_mi_s0_e32 = a0*s1*inv_mi_s0_e32;

	complex<double> a2 = (m_AuxIntPar.half_sqrt_pi*inv_mi_s1_e32)*(((m_AuxIntPar.two_sqrt_pi*m_AuxIntPar.pxz)*a0_s1_inv_mi_s0_e32) + a1*sx1z); //OC26042017 (uncommented)
		//complex<double> a2 = 0;

	complex<double> pi_d_sqrt_mi_s0_d_sqrt_mi_s1 = m_AuxIntPar.pi/(sqrt_mi_s0*sqrt_mi_s1);
	complex<double> a2z = pi_d_sqrt_mi_s0_d_sqrt_mi_s1*vz;
	complex<double> half_sqrt_pi_inv_mi_s2_e32 = m_AuxIntPar.half_sqrt_pi*inv_mi_s2_e32;
	
	complex<double> a3 = half_sqrt_pi_inv_mi_s2_e32*(((-m_AuxIntPar.sqrt_pi)*s2*inv_mi_s1_e32)*(((m_AuxIntPar.two_sqrt_pi*m_AuxIntPar.pxz1)*a0_s1_inv_mi_s0_e32) + (a1*sx1z1)) + (a2*szz1a)); //OC26042017 (uncommented)
		//complex<double> a3 = 0;

	complex<double> a3z = half_sqrt_pi_inv_mi_s2_e32*((a2z*szz1a) - (2.*s2*vz1*pi_d_sqrt_mi_s0_d_sqrt_mi_s1));
	complex<double> inv_mi_s2_e32_inv_mi_s3_e32 = inv_mi_s2_e32*inv_mi_s3_e32;
    
	complex<double> a4 = (0.5*inv_mi_s0_e32*inv_mi_s1_e32*inv_mi_s2_e32_inv_mi_s3_e32)*(((2*m_AuxIntPar.piE2*m_AuxIntPar.pxs)*a0*s1*s2*s3) + (m_AuxIntPar.sqrt_pi*mi_s0_e32)*((m_AuxIntPar.pi*a1*s2*s3*sx1s) + (sqrt_mi_s1*s1)*((a3*sqrt_mi_s2*s2*sz1s) + (a2*m_AuxIntPar.sqrt_pi*s3*szs)))); //OC26042017 (uncommented)
		//complex<double> a4 = 0;
	
	complex<double> a4z = (m_AuxIntPar.half_sqrt_pi*a3z*sz1s*inv_mi_s3_e32) - (0.5*inv_mi_s2_e32_inv_mi_s3_e32*m_AuxIntPar.pi*a2z*s3*szs); //OC26042017 (uncommented)
		//complex<double> a4z = 0;
	
	complex<double> half_sqrt_pi_inv_mi_s4_e32 = m_AuxIntPar.half_sqrt_pi*inv_mi_s4_e32;
	complex<double> mi_sqrt_pi_s4_inv_mi_s3_e32 = -m_AuxIntPar.sqrt_pi*s4*inv_mi_s3_e32;
    complex<double> mi_sqrt_pi_s3_inv_mi_s2_e32 = -m_AuxIntPar.sqrt_pi*s3*inv_mi_s2_e32;
    complex<double> mi_sqrt_pi_s2_inv_mi_s1_e32 = -m_AuxIntPar.sqrt_pi*s2*inv_mi_s1_e32;
	complex<double> mi_sqrt_pi_s1_inv_mi_s0_e32 = -m_AuxIntPar.sqrt_pi*s1*inv_mi_s0_e32;
	
	complex<double> a5 = half_sqrt_pi_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*((-m_AuxIntPar.sqrt_pi*s1*inv_mi_s0_e32)*((-2.*hg*s0) + (a0*sxg)) + (a1*sx1g)) + (a2*szga)) + (a3*sz1g)) + (a4*ssg)); //OC26042017 (uncommented)
		//complex<double> a5 = half_sqrt_pi_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*((-m_AuxIntPar.sqrt_pi*s1*inv_mi_s0_e32)*((-2.*hg*s0) + (a0*sxg)) + (a1*sx1g)))));
	
	complex<double> a5z = half_sqrt_pi_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(-m_AuxIntPar.two_sqrt_pi*s1*vg/sqrt_mi_s0) + (a2z*szga)) + (a3z*sz1g)) + (a4z*ssg)); //OC26042017 (uncommented)
		//complex<double> a5z = half_sqrt_pi_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(-m_AuxIntPar.two_sqrt_pi*s1*vg/sqrt_mi_s0) + (a2z*szga)) + (a3z*sz1g)));

    complex<double> w1(m_AuxIntPar.constW1, phi_x);
	complex<double> w1_d_s0 = w1/s0;
	complex<double> w2(m_AuxIntPar.constW2, phi_x1);
	w2 -= 0.5*w1_d_s0*sxx1;
	complex<double> half_w2_d_s1 = 0.5*w2/s1;
	complex<double> w3a(m_AuxIntPar.constW3a, phi_z);
    
	w3a += (m_AuxIntPar.pxz*w1_d_s0) - (sx1z*half_w2_d_s1); //OC26042017 (uncommented)

	complex<double> half_w3a_d_s2 = 0.5*w3a/s2;
	complex<double> w3b(m_AuxIntPar.constW3b, phi_z1);

	w3b += (m_AuxIntPar.pxz1*w1_d_s0) - (sx1z1*half_w2_d_s1); //OC26042017 (uncommented)

	complex<double> w4a = w3b - (half_w3a_d_s2*szz1a);
	complex<double> half_w4a_d_s3 = 0.5*w4a/s3;
	complex<double> w4b(m_AuxIntPar.constW4b, phi_s);

	w4b += (m_AuxIntPar.pxs*w1_d_s0) - (sx1s*half_w2_d_s1) - (szs*half_w3a_d_s2) - (sz1s*half_w4a_d_s3); //OC26042017 (uncommented)
	
	complex<double> half_w4b_d_s4 = 0.5*w4b/s4;
	complex<double> w5(m_AuxIntPar.constW5, phi_g);

	w5 -= (0.5*sxg*w1_d_s0) + (sx1g*half_w2_d_s1) + (szga*half_w3a_d_s2) + (sz1g*half_w4a_d_s3) + (ssg*half_w4b_d_s4); //OC26042017 (uncommented)
		//w5 -= (0.5*sxg*w1_d_s0) + (sx1g*half_w2_d_s1) + (szga*half_w3a_d_s2) + (sz1g*half_w4a_d_s3);

	complex<double> half_sqrt_pi_inv_mi_s5_e32 = m_AuxIntPar.half_sqrt_pi*inv_mi_s5_e32;
	complex<double> mi_sqrt_pi_s5_inv_mi_s4_e32 = -m_AuxIntPar.sqrt_pi*s5*inv_mi_s4_e32;
	
	ampX = half_sqrt_pi_inv_mi_s5_e32*(mi_sqrt_pi_s5_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(mi_sqrt_pi_s1_inv_mi_s0_e32*((-2.*h0*s0) + (a0*w1)) + (a1*w2)) + (a2*w3a)) + (a3*w4a)) + (a4*w4b)) + (a5*w5)); //OC26042017 (uncommented)
		//ampX = half_sqrt_pi_inv_mi_s5_e32*(mi_sqrt_pi_s5_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(mi_sqrt_pi_s1_inv_mi_s0_e32*((-2.*h0*s0) + (a0*w1)) + (a1*w2))))) + (a5*w5));
	
	ampZ = half_sqrt_pi_inv_mi_s5_e32*(mi_sqrt_pi_s5_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(mi_sqrt_pi_s1_inv_mi_s0_e32*(-2.*s0*v0)) + (a2z*w3a)) + (a3z*w4a)) + (a4z*w4b)) + (a5z*w5)); //OC26042017 (uncommented)
		//ampZ = half_sqrt_pi_inv_mi_s5_e32*(mi_sqrt_pi_s5_inv_mi_s4_e32*(mi_sqrt_pi_s4_inv_mi_s3_e32*(mi_sqrt_pi_s3_inv_mi_s2_e32*(mi_sqrt_pi_s2_inv_mi_s1_e32*(mi_sqrt_pi_s1_inv_mi_s0_e32*(-2.*s0*v0)) + (a2z*w3a)) + (a3z*w4a))) + (a5z*w5));
	
	complex<double> q6(m_AuxIntPar.constQ6, phi_0);

	complex<double> w4b_half_w4b_d_s4 = w4b*half_w4b_d_s4;

	q6 -= 0.5*((0.5*w1*w1_d_s0) + (w2*half_w2_d_s1) + (w3a*half_w3a_d_s2) + (w4a*half_w4a_d_s3) + w4b_half_w4b_d_s4 + (0.5*w5*w5/s5));

	complex<double> extraArg(0, k*exzy.y + m_AuxIntPar.half_pi);
    arg = q6 + extraArg;

	double arg_re = arg.real(), arg_im = arg.imag();
	double exp_arg_re = exp(arg_re), cos_arg_im = cos(arg_im), sin_arg_im = sin(arg_im);

	//complex<double> expQ6 = exp(arg);
	complex<double> expQ6(exp_arg_re*cos_arg_im, exp_arg_re*sin_arg_im);

	complex<double> FuncX = ampX*expQ6, FuncZ = ampZ*expQ6;

    Ew.EwX_Re = FuncX.real(); Ew.EwX_Im = FuncX.imag();
    Ew.EwZ_Re = FuncZ.real(); Ew.EwZ_Im = FuncZ.imag();
}

//*************************************************************************

//void srTCSR::createAuxWfrSmp()
//{
//	m_pWfrSmpAux = new srTWfrSmp(m_Wfr.yStart, m_Wfr.xStart, m_Wfr.xStart + m_Wfr.xStep*(m_Wfr.nx - 1), m_Wfr.nx, m_Wfr.zStart, m_Wfr.zStart + m_Wfr.zStep*(m_Wfr.nz - 1), m_Wfr.nz, m_Wfr.eStart, m_Wfr.eStart + m_Wfr.eStep*(m_Wfr.ne - 1), m_Wfr.ne, "EV");
//}

//*************************************************************************

void srTCSR::radIntegrationMonteCarlo(srTEXZY& exzy, srTEFourier& Ew)
{//dedicated for test
	Ew.EwX_Re = Ew.EwX_Im = Ew.EwZ_Re = Ew.EwZ_Im = 0.;
	srTEbmDat& e_beam = m_pTrjDatAux->EbmDat;

    double k = exzy.e*m_AuxIntPar.k_d_e;
	char rand_mode = 1; //use LPTau
	double point6d[6];
	for(int i=0; i<m_PrecParams.m_MC_NumMacroPart; i++)
	{
		m_pWfrSmpAux->yStart = m_pWfrSmpAux->yEnd = exzy.y;
		m_pWfrSmpAux->zStart = m_pWfrSmpAux->zEnd = exzy.z;
		m_pWfrSmpAux->xStart = m_pWfrSmpAux->xEnd = exzy.x;
		m_pWfrSmpAux->LambStart = m_pWfrSmpAux->LambEnd = exzy.e;
		m_pWfrAux->yStart = exzy.y;
		m_pWfrAux->zStart = exzy.z;
		m_pWfrAux->xStart = exzy.x;
		m_pWfrAux->eStart = exzy.e;

		float *pEx = m_pWfrAux->pBaseRadX, *pEz = m_pWfrAux->pBaseRadZ;
		*pEx = 0; *(pEx + 1) = 0; *pEz = 0; *(pEz + 1) = 0;
		
		m_gmRand.NextRandGauss6D(m_xcArr, m_sigArr, point6d, rand_mode);
		e_beam.x0 = point6d[0];
		e_beam.dxds0 = point6d[1];
		e_beam.z0 = point6d[2];
		e_beam.dzds0 = point6d[3];
		e_beam.SetNewEnergy(point6d[4]);
		e_beam.sc = point6d[5];

		m_pTrjDatAux->ComputeInterpolatingStructure();
		m_pRadInt->ComputeElectricFieldFreqDomain(m_pTrjDatAux, m_pWfrSmpAux, m_pParPrecElecFldSingle, m_pWfrAux, 0);

		double ksc = -k*e_beam.sc;
        double cos_ksc = cos(ksc), sin_ksc = sin(ksc);

		double newEwX_Re = (*pEx)*cos_ksc - (*(pEx + 1))*sin_ksc;
		double newEwX_Im = (*pEx)*sin_ksc + (*(pEx + 1))*cos_ksc;
		double newEwZ_Re = (*pEz)*cos_ksc - (*(pEz + 1))*sin_ksc;
		double newEwZ_Im = (*pEz)*sin_ksc + (*(pEz + 1))*cos_ksc;

		m_pRadInt->DeallocateMemForRadDistr();

		double inv_i_p_1 = 1./(double(i + 1));
		Ew.EwX_Re = (i*Ew.EwX_Re + newEwX_Re)*inv_i_p_1;
		Ew.EwX_Im = (i*Ew.EwX_Im + newEwX_Im)*inv_i_p_1;
		Ew.EwZ_Re = (i*Ew.EwZ_Re + newEwZ_Re)*inv_i_p_1;
		Ew.EwZ_Im = (i*Ew.EwZ_Im + newEwZ_Im)*inv_i_p_1;

			//if((i==200) || (i==400) || (i==600) || (i==800) || (i==1000) || (i==1200) || (i==1400) || (i==1600) || (i==1800) || (i==2000))
			//{
			//	int aha = 1;
			//}
	}
	Ew *= sqrt(e_beam.Neb);
}

//*************************************************************************
