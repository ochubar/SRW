/************************************************************************//**
 * File: sroptcryst.cpp
 * Description: Optical element: Perfect Crystal
 * Project: Synchrotron Radiation Workshop
 * First release: October 2014
 *
 * Copyright (C) Brookhaven National Laboratory, Upton, NY, USA
 * Portions Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#include "sroptcryst.h"

//*************************************************************************

int srTOptCryst::WfrInterpolOnRegGrid(srTSRWRadStructAccessData* pWfr, srTOptCrystMeshTrf* pMeshTrf)
{
	if((pWfr == 0) || (pMeshTrf == 0)) return FAILED_INTERPOL_ELEC_FLD;

	char quadWfrTermCanBeTreated = WaveFrontTermCanBeTreated(*pWfr, false);
	bool isCoordRepres = (pWfr->Pres == 0);
	bool isFreqRepres = (pWfr->PresT == 0);

	//long PerE = 2;
	//long PerX = PerE << 1;
	//long PerZ = PerX*(pWfr->nx);

	long long PerE = 2;
	long long PerX = PerE << 1;
	long long PerZ = PerX*(pWfr->nx);

	//float *t_ExRes = pWfr->pBaseRadX;
	//float *t_EzRes = pWfr->pBaseRadZ;

	//long nTot = PerZ*(pWfr->nz);
	//long nTot_mi_1 = nTot - 1;
	//long nx_mi_1 = pWfr->nx - 1;
	//long nz_mi_1 = pWfr->nz - 1;
	//double f0m1, fm10, f00, f10, f01, f11, a10, a01, a11, a20, a02;

	//long ix0=-1, iz0=-1;
	//double phEn, x, z = pWfr->zStart;

	//const double dMax = 1.E+20;
	//double dx, dz, dTest, dTest0, dTestPrev;

	for(long iz=0; iz<pWfr->nz; iz++)
	{
		//ddddddddddddddddddddddddddddddddddddddddddddddd
/**
		x = pWfr->xStart;
		for(long ix=0; ix<pWfr->nx; ix++)
		{
			bool pointIsInsideNonZeroSquare = (x >= xMin) && (x <= xMax) && (z >= zMin) && (z <= zMax);

			phEn = pWfr->eStart;
			for(long ie=0; ie<pWfr->ne; ie++)
			{
				long two_ie = 2*ie;

				if(pointIsInsideNonZeroSquare)
				{//perform interpolation on irregular mesh (bilinear or based on 12 points)
					//find indexes of the relevant point for the interpolation

					//OCTEST
					//if((fabs(pWfr->RobsX + 0.99) < 0.1) && (fabs(pWfr->RobsZ + 0.99) < 0.1))
					//{
					//	//if((fabs(rx_m10) > dMax) || (fabs(rx_10) > dMax) || (fabs(rz_0m1) > dMax) || (fabs(rz_01) > dMax))
					//	//if((iz == 1848) && (ix == 442))// && (fabs(x + 0.0047703) < 0.00001))
					//	if((iz == 748) && (ix == 350))// && (fabs(x + 0.0047703) < 0.00001))
					//	{ 
					//		int aha = 1;
					//	}
					//}

					if(ix0 < 0) ix0 = ix;
					if(iz0 < 0) iz0 = iz;

					bool pointFound = false, candPointFound = false;
					//bool isLeftBordX = false, isRightBordX = false;
					//bool isLeftBordZ = false, isRightBordZ = false;

					long ix0pr = -1, iz0pr = -1;
					while((ix0 != ix0pr) && (iz0 != iz0pr)) 
					{//This while loop is required for a "tilted/rotated" mesh  (to check how ir works!)
						ix0pr = ix0; iz0pr = iz0;

						long iz0_PerZ_p_2_ie = iz0*PerZ + two_ie;

						dTestPrev = 1.E+23;
						long ofst = ix0*PerX + iz0_PerZ_p_2_ie;
						if (ofst < nTot_mi_1)
						{
							dx = x - arRayTrCoord[ofst];
							dz = z - arRayTrCoord[ofst + 1];
							dTestPrev = sqrt(dx*dx + dz*dz);
						}
						dTest0 = dTestPrev;

						pointFound = false;
						candPointFound = false;
						//isLeftBordX = false; isRightBordX = false;

						long ix0orig = ix0;
						//for(int iix = ix0 - 1; iix >= 0; iix--)
						for(int iix = ix0orig - 1; iix >= 0; iix--) //OC200414
						{
							ofst = iix*PerX + iz0_PerZ_p_2_ie;
							if(ofst < nTot_mi_1)
							{
								dx = x - arRayTrCoord[ofst];
								dz = z - arRayTrCoord[ofst + 1];
								dTest = sqrt(dx*dx + dz*dz);
								//if((dTest > dMax) && (!candPointFound)) continue;
								if(dTest > dMax)
								{
									if(dTestPrev < dMax) break;
									if(!candPointFound) continue;
								}
								if(dTest < dTestPrev)
								{
									ix0 = iix; dTestPrev = dTest; candPointFound = true;
								}
								else
								{
									if(candPointFound)
									{
										pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isLeftBordX = true;
									}
									break;
								}
							}
						}
						if(!pointFound)
						{
							//dTestPrev = 1.E+23;
							dTestPrev = dTest0;
							candPointFound = false;
							//for(int iix = ix0 + 1; iix < pWfr->nx; iix++)
							for(int iix = ix0orig + 1; iix < pWfr->nx; iix++) //OC200414
							{
								ofst = iix*PerX + iz0_PerZ_p_2_ie;
								if(ofst < nTot_mi_1)
								{
									dx = x - arRayTrCoord[ofst];
									dz = z - arRayTrCoord[ofst + 1];
									dTest = sqrt(dx*dx + dz*dz);
									//if((dTest > dMax) && (!candPointFound)) continue;
									if(dTest > dMax)
									{
										if(dTestPrev < dMax) break;
										if(!candPointFound) continue;
									}
									if(dTest < dTestPrev)
									{
										ix0 = iix; dTestPrev = dTest; candPointFound = true;
									}
									else
									{
										//if(candPointFound) pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isRightBordX = true;
										break;
									}
								}
							}
						}

						long ix0_PerX_p_2_ie = ix0*PerX + two_ie;

						dTestPrev = 1.E+23;
						ofst = iz0*PerZ + ix0_PerX_p_2_ie;
						if (ofst < nTot_mi_1)
						{
							dx = x - arRayTrCoord[ofst];
							dz = z - arRayTrCoord[ofst + 1];
							dTestPrev = sqrt(dx*dx + dz*dz);
						}
						dTest0 = dTestPrev;

						pointFound = false;
						candPointFound = false;
						//isLeftBordZ = false; isRightBordZ = false;

						long iz0orig = iz0;
						//for(int iiz = iz0 - 1; iiz >= 0; iiz--)
						for(int iiz = iz0orig - 1; iiz >= 0; iiz--) //OC200414
						{
							ofst = iiz*PerZ + ix0_PerX_p_2_ie;
							if(ofst < nTot_mi_1)
							{
								dx = x - arRayTrCoord[ofst];
								dz = z - arRayTrCoord[ofst + 1];
								dTest = sqrt(dx*dx + dz*dz);
								//if((dTest > dMax) && (!candPointFound)) continue;
								if(dTest > dMax)
								{
									if(dTestPrev < dMax) break;
									if(!candPointFound) continue;
								}
								if(dTest < dTestPrev)
								{
									iz0 = iiz; dTestPrev = dTest; candPointFound = true;
								}
								else
								{
									if(candPointFound)
									{
										pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isLeftBordZ = true;
									}
									//OCTEST
									//if((dTest > dMax) || (dTestPrev > dMax))
									//{
									//	int aha = 1;
									//}
									break;
								}
							}
						}
						if(!pointFound)
						{
							//dTestPrev = 1.E+23;
							dTestPrev = dTest0;
							candPointFound = false;
							//for(int iiz = iz0 + 1; iiz < pWfr->nz; iiz++)
							for(int iiz = iz0orig + 1; iiz < pWfr->nz; iiz++) //OC200414
							{
								ofst = iiz*PerZ + ix0_PerX_p_2_ie;
								if(ofst < nTot_mi_1)
								{
									dx = x - arRayTrCoord[ofst];
									dz = z - arRayTrCoord[ofst + 1];
									dTest = sqrt(dx*dx + dz*dz);
									//if((dTest > dMax) && (!candPointFound)) continue;
									if(dTest > dMax)
									{
										if(dTestPrev < dMax) 
										{
											//if(dx < 0) isRightBordZ = true;
											//else if(dx > 0) isLeftBordZ = true;
											break;
										}
										if(!candPointFound) continue;
									}
									if(dTest < dTestPrev)
									{
										iz0 = iiz; dTestPrev = dTest; candPointFound = true;
									}
									else
									{
										//if(candPointFound) pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) 
										//{
										//	if(dx < 0) isRightBordZ = true;
										//	else if(dx > 0) isLeftBordZ = true;
										//}

										//OCTEST
										//if((dTest > dMax) || (dTestPrev > dMax))
										//{
										//	int aha = 1;
										//}
										break;
									}
								}
							}
						}
					}
					//calculate indexes of other points and interpolate 
					//2 cases are considered: "bi-linear" 2D interpolation based on 4 points and "bi-quadratic" 2D interpolation based on 5 points (mesh can be irregular)
					const double relTolEqualStep = 1.e-04; //to tune
					
					if(ix0 < 0) ix0 = 0;
					else if(ix0 >= nx_mi_1) ix0 = nx_mi_1 - 1;
					if(iz0 < 0) iz0 = 0;
					else if(iz0 >= nz_mi_1) iz0 = nz_mi_1 - 1;

					long ofst_00 = ix0*PerX + iz0*PerZ + two_ie;
					long ofst_m10 = ofst_00 - PerX;
					long ofst_10 = ofst_00 + PerX;
					long ofst_0m1 = ofst_00 - PerZ;
					long ofst_01 = ofst_00 + PerZ;

					if(ix0 == 0) ofst_m10 = ofst_00;
					if(ix0 == nx_mi_1) ofst_10 = ofst_00;
					if(iz0 == 0) ofst_0m1 = ofst_00;
					if(iz0 == nz_mi_1) ofst_10 = ofst_00;

					//if((ix0 == 0) || isLeftBordX) ofst_m10 = ofst_00;
					//if((ix0 == nx_mi_1) || isRightBordX) ofst_10 = ofst_00;
					//if((iz0 == 0) || isLeftBordZ) ofst_0m1 = ofst_00;
					//if((iz0 == nz_mi_1) || isRightBordZ) ofst_10 = ofst_00;

					long ofst_00_p_1 = ofst_00 + 1;
					long ofst_m10_p_1 = ofst_m10 + 1;
					long ofst_10_p_1 = ofst_10 + 1;
					long ofst_0m1_p_1 = ofst_0m1 + 1;
					long ofst_01_p_1 = ofst_01 + 1;

					double x_00 = arRayTrCoord[ofst_00], z_00 = arRayTrCoord[ofst_00_p_1];
					double x_m10 = arRayTrCoord[ofst_m10], z_m10 = arRayTrCoord[ofst_m10_p_1];
					double x_10 = arRayTrCoord[ofst_10], z_10 = arRayTrCoord[ofst_10_p_1];
					double x_0m1 = arRayTrCoord[ofst_0m1], z_0m1 = arRayTrCoord[ofst_0m1_p_1];
					double x_01 = arRayTrCoord[ofst_01], z_01 = arRayTrCoord[ofst_01_p_1];

					double rx_m10 = x_m10 - x_00, rz_m10 = z_m10 - z_00;
					double rx_10 = x_10 - x_00, rz_10 = z_10 - z_00;
					double rx_0m1 = x_0m1 - x_00, rz_0m1 = z_0m1 - z_00;
					double rx_01 = x_01 - x_00, rz_01 = z_01 - z_00;
					double dx_00 = x - x_00, dz_00 = z - z_00;

					//OCTEST
					//if(fabs(pWfr->RobsZ + 0.37) < 0.1)
					//{
					//	//if((fabs(rx_m10) > dMax) || (fabs(rx_10) > dMax) || (fabs(rz_0m1) > dMax) || (fabs(rz_01) > dMax))
					//	if((fabs(z + 0.0054932) < 0.00001)) //&& (fabs(x + 0.0047703) < 0.00001))
					//	{
					//		int aha = 1;
					//	}
					//}

					bool rx_m10_isNotOK = ((fabs(rx_m10) > dMax) || (rx_m10 == 0));
					bool rx_10_isNotOK = ((fabs(rx_10) > dMax) || (rx_10 == 0));
					if(rx_m10_isNotOK && rx_10_isNotOK) goto SetFieldToZero;

					bool rz_0m1_isNotOK = ((fabs(rz_0m1) > dMax) || (rz_0m1 == 0));
					bool rz_01_isNotOK = ((fabs(rz_01) > dMax) || (rz_01 == 0));
					if(rz_0m1_isNotOK && rz_01_isNotOK) goto SetFieldToZero;

					if(rx_m10_isNotOK) 
					{
						rx_m10 = -rx_10;
						ofst_m10 = ofst_00;
					}
					else if(rx_10_isNotOK) 
					{
						rx_10 = -rx_m10;
						ofst_10 = ofst_00;
					}

					bool rx_0m1_isNotOK = (fabs(rx_0m1) > dMax);
					bool rx_01_isNotOK = (fabs(rx_01) > dMax);
					if(rx_0m1_isNotOK) 
					{
						rx_0m1 = 0.; //??
						ofst_0m1 = ofst_00;
					}
					if(rx_01_isNotOK) 
					{
						rx_01 = 0.;
						ofst_01 = ofst_00;
					}

					if(rz_0m1_isNotOK) 
					{
						rz_0m1 = -rz_01;
						ofst_0m1 = ofst_00;
					}
					else if(rz_01_isNotOK) 
					{
						rz_01 = -rz_0m1;
						ofst_01 = ofst_00;
					}

					bool rz_m10_isNotOK = (fabs(rz_m10) > dMax);
					bool rz_10_isNotOK = (fabs(rz_10) > dMax);
					if(rz_m10_isNotOK) 
					{
						rz_m10 = 0.;
						ofst_m10 = ofst_00;
					}
					if(rz_10_isNotOK) 
					{
						rz_10 = 0.;
						ofst_10 = ofst_00;
					}

					const double maxRelArg = 1.5;
					//OC200414
					//if(rx_m10_isNotOK || rx_10_isNotOK || rx_0m1_isNotOK || rx_01_isNotOK)
					//{
					//	double twoRx = ::fabs(rx_m10) + ::fabs(rx_10);
					//	if(::fabs(dx_00/twoRx) > maxRelArg) goto SetFieldToZero;
					//}
					//if(rz_0m1_isNotOK || rz_01_isNotOK || rz_m10_isNotOK || rz_10_isNotOK)
					//{
					//	double twoRz = ::fabs(rz_0m1) + ::fabs(rz_01);
					//	if(::fabs(dz_00/twoRz) > maxRelArg) goto SetFieldToZero;
					//}

					double d_rx_m10 = dx_00 - rx_m10;
					double d_rx_10 = dx_00 - rx_10;
					if(d_rx_m10*d_rx_10 > 0.)
					{
						if((::fabs(rx_m10)*maxRelArg < ::fabs(d_rx_m10)) || (::fabs(rx_10)*maxRelArg < ::fabs(d_rx_10))) goto SetFieldToZero;
					}
					double d_rz_0m1 = dz_00 - rz_0m1;
					double d_rz_01 = dz_00 - rz_01;
					if(d_rz_0m1*d_rz_01 > 0.)
					{
						if((::fabs(rz_0m1)*maxRelArg < ::fabs(d_rz_0m1)) || (::fabs(rz_01)*maxRelArg < ::fabs(d_rz_01))) goto SetFieldToZero;
					}

					if(m_wfrInterpMode == 1)
					{//bi-linear, based on 4 points
						double sp_m10 = rx_m10*dx_00 + rz_m10*dz_00;
						double sp_10 = rx_10*dx_00 + rz_10*dz_00;
						double sp_0m1 = rx_0m1*dx_00 + rz_0m1*dz_00;
						double sp_01 = rx_01*dx_00 + rz_01*dz_00;

						bool initPointMoved = false;
						if((sp_m10 > 0) && (sp_10 < 0))
						{
							if(ix0 > 0) { ix0--; initPointMoved = true;}
						}
						if((sp_0m1 > 0) && (sp_01 < 0))
						{
							if(iz0 > 0) { iz0--; initPointMoved = true;}
						}

						long ofst0 = initPointMoved? (ix0*PerX + iz0*PerZ + two_ie) : ofst_00;

						long ofst0_p_PerX = ofst0 + PerX;
						long ofst0_p_PerZ = ofst0 + PerZ;
						long ofst0_p_PerX_p_PerZ = ofst0_p_PerZ + PerX;
						double x00 = arRayTrCoord[ofst0], x10 = arRayTrCoord[ofst0_p_PerX], x01 = arRayTrCoord[ofst0_p_PerZ], x11 = arRayTrCoord[ofst0_p_PerX_p_PerZ];
						long ofst0_p_1 = ofst0 + 1;
						long ofst0_p_PerX_p_1 = ofst0_p_PerX + 1;
						long ofst0_p_PerZ_p_1 = ofst0_p_PerZ + 1;
						long ofst0_p_PerX_p_PerZ_p_1 = ofst0_p_PerX_p_PerZ + 1;
						double z00 = arRayTrCoord[ofst0_p_1], z10 = arRayTrCoord[ofst0_p_PerX_p_1], z01 = arRayTrCoord[ofst0_p_PerZ_p_1], z11 = arRayTrCoord[ofst0_p_PerX_p_PerZ_p_1];

						double rX = x - x00, rZ = z - z00;
						double rX10 = x10 - x00, rZ10 = z10 - z00;
						double rX01 = x01 - x00, rZ01 = z01 - z00;
						double rX11 = x11 - x00, rZ11 = z11 - z00;

						double absTolX = relTolEqualStep*fabs(rX10);
						double absTolZ = relTolEqualStep*fabs(rZ01);
						bool isRecX = (fabs(rX01) < absTolX) && (fabs(x11 - x10) < absTolX);
						bool isRecZ = (fabs(rZ10) < absTolZ) && (fabs(z11 - z01) < absTolZ);
						if(isRecX && isRecZ)
						{//regular rectangular mesh
							double xt = rX/rX10, zt = rZ/rZ01;
							if(arEX != 0)
							{
								f00 = arEX[ofst0]; f10 = arEX[ofst0_p_PerX]; f01 = arEX[ofst0_p_PerZ]; f11 = arEX[ofst0_p_PerX_p_PerZ];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*t_ExRes = (float)(xt*(a10 + a11*zt) + a01*zt + f00);

								f00 = arEX[ofst0_p_1]; f10 = arEX[ofst0_p_PerX_p_1]; f01 = arEX[ofst0_p_PerZ_p_1]; f11 = arEX[ofst0_p_PerX_p_PerZ_p_1];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*(t_ExRes + 1) = (float)(xt*(a10 + a11*zt) + a01*zt + f00);
							}
							if(arEZ != 0)
							{
								f00 = arEZ[ofst0]; f10 = arEZ[ofst0_p_PerX]; f01 = arEZ[ofst0_p_PerZ]; f11 = arEZ[ofst0_p_PerX_p_PerZ];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*t_EzRes = (float)(xt*(a10 + a11*zt) + a01*zt + f00);

								f00 = arEZ[ofst0_p_1]; f10 = arEZ[ofst0_p_PerX_p_1]; f01 = arEZ[ofst0_p_PerZ_p_1]; f11 = arEZ[ofst0_p_PerX_p_PerZ_p_1];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*(t_EzRes + 1) = (float)(xt*(a10 + a11*zt) + a01*zt + f00);
							}
						}
						else
						{//irregular mesh (general case)
							double arXZ[] = {rX10, rZ10, rX01, rZ01, rX11, rZ11};
							if(arEX != 0)
							{
								double arER[] = {arEX[ofst0], arEX[ofst0_p_PerX], arEX[ofst0_p_PerZ], arEX[ofst0_p_PerX_p_PerZ]};
								*t_ExRes = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arER);
								double arEI[] = {arEX[ofst0_p_1], arEX[ofst0_p_PerX_p_1], arEX[ofst0_p_PerZ_p_1], arEX[ofst0_p_PerX_p_PerZ_p_1]};
								*(t_ExRes + 1) = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arEI);
							}
							if(arEZ != 0)
							{
								double arER[] = {arEZ[ofst0], arEZ[ofst0_p_PerX], arEZ[ofst0_p_PerZ], arEZ[ofst0_p_PerX_p_PerZ]};
								*t_EzRes = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arER);
								double arEI[] = {arEZ[ofst0_p_1], arEZ[ofst0_p_PerX_p_1], arEZ[ofst0_p_PerZ_p_1], arEZ[ofst0_p_PerX_p_PerZ_p_1]};
								*(t_EzRes + 1) = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arEI);
							}
						}
					}
					else if(m_wfrInterpMode == 2)
					{//bi-quadratic, based on 5 points
						double absTolX = relTolEqualStep*fabs(rx_10);
						double absTolZ = relTolEqualStep*fabs(rz_01);
						bool isRecX = (fabs(rx_0m1) < absTolX) && (fabs(rx_01) < absTolX);
						bool isRecZ = (fabs(rz_m10) < absTolZ) && (fabs(rz_10) < absTolX);
						if(isRecX && isRecZ)
						{
							bool isEquidistX = (fabs(rx_m10 + rx_10) < absTolX);
							bool isEquidistZ = (fabs(rz_0m1 + rz_01) < absTolZ);
							if(isEquidistX && isEquidistZ)
							{//regular rectangular mesh
								double xt = dx_00/rx_10, zt = dz_00/rz_01;
								if(arEX != 0)
								{
									f0m1 = arEX[ofst_0m1]; fm10 = arEX[ofst_m10]; f00 = arEX[ofst_00]; f10 = arEX[ofst_10]; f01 = arEX[ofst_01];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*t_ExRes = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);

									f0m1 = arEX[ofst_0m1_p_1]; fm10 = arEX[ofst_m10_p_1]; f00 = arEX[ofst_00_p_1]; f10 = arEX[ofst_10_p_1]; f01 = arEX[ofst_01_p_1];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*(t_ExRes + 1) = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);
								}
								if(arEZ != 0)
								{
									f0m1 = arEZ[ofst_0m1]; fm10 = arEZ[ofst_m10]; f00 = arEZ[ofst_00]; f10 = arEZ[ofst_10]; f01 = arEZ[ofst_01];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*t_EzRes = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);

									f0m1 = arEZ[ofst_0m1_p_1]; fm10 = arEZ[ofst_m10_p_1]; f00 = arEZ[ofst_00_p_1]; f10 = arEZ[ofst_10_p_1]; f01 = arEZ[ofst_01_p_1];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*(t_EzRes + 1) = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);
								}
							}
							else
							{//variable-step rectangular mesh
								double DX = 1./((rx_10 - rx_m10)*rx_10*rx_m10), DZ = 1./((rz_01 - rz_0m1)*rz_01*rz_0m1);
								if(arEX != 0)
								{
									f0m1 = arEX[ofst_0m1]; fm10 = arEX[ofst_m10]; f00 = arEX[ofst_00]; f10 = arEX[ofst_10]; f01 = arEX[ofst_01];
									double f00_mi_fm10_x1 = (f00 - fm10)*rx_10, f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									double f00_mi_f10_xm1 = (f00 - f10)*rx_m10, f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*t_ExRes = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);

									f0m1 = arEX[ofst_0m1_p_1]; fm10 = arEX[ofst_m10_p_1]; f00 = arEX[ofst_00_p_1]; f10 = arEX[ofst_10_p_1]; f01 = arEX[ofst_01_p_1];
									f00_mi_fm10_x1 = (f00 - fm10)*rx_10; f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									f00_mi_f10_xm1 = (f00 - f10)*rx_m10; f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*(t_ExRes + 1) = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);
								}
								if(arEZ != 0)
								{
									f0m1 = arEZ[ofst_0m1]; fm10 = arEZ[ofst_m10]; f00 = arEZ[ofst_00]; f10 = arEZ[ofst_10]; f01 = arEZ[ofst_01];
									double f00_mi_fm10_x1 = (f00 - fm10)*rx_10, f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									double f00_mi_f10_xm1 = (f00 - f10)*rx_m10, f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*t_EzRes = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);

									f0m1 = arEZ[ofst_0m1_p_1]; fm10 = arEZ[ofst_m10_p_1]; f00 = arEZ[ofst_00_p_1]; f10 = arEZ[ofst_10_p_1]; f01 = arEZ[ofst_01_p_1];
									f00_mi_fm10_x1 = (f00 - fm10)*rx_10; f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									f00_mi_f10_xm1 = (f00 - f10)*rx_m10; f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*(t_EzRes + 1) = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);
								}
							}
						}
						else
						{//irregular mesh (general case)
							double arXZ[] = {rx_0m1, rz_0m1, rx_m10, rz_m10, rx_10, rz_10, rx_01, rz_01};
							if(arEX != 0)
							{
								double arER[] = {arEX[ofst_0m1], arEX[ofst_m10], arEX[ofst_00], arEX[ofst_10], arEX[ofst_01]};
								*t_ExRes = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arER);
								double arEI[] = {arEX[ofst_0m1_p_1], arEX[ofst_m10_p_1], arEX[ofst_00_p_1], arEX[ofst_10_p_1], arEX[ofst_01_p_1]};
								*(t_ExRes + 1) = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arEI);
							}
							if(arEZ != 0)
							{
								double arER[] = {arEZ[ofst_0m1], arEZ[ofst_m10], arEZ[ofst_00], arEZ[ofst_10], arEZ[ofst_01]};
								*t_EzRes = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arER);
								double arEI[] = {arEZ[ofst_0m1_p_1], arEZ[ofst_m10_p_1], arEZ[ofst_00_p_1], arEZ[ofst_10_p_1], arEZ[ofst_01_p_1]};
								*(t_EzRes + 1) = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arEI);
							}
						}
					}
				}
				else
				{
					SetFieldToZero:
					*t_ExRes = 0.; *(t_ExRes+1) = 0.;
					*t_EzRes = 0.; *(t_EzRes+1) = 0.;

					ix0 = -1; iz0 = -1; //OC200414 (forces to restart search)
				}

				//test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//long ofstAux = iz*PerZ + ix*PerX + ie*2;
				//*t_ExRes = arEX[ofstAux]; *(t_ExRes+1) = arEX[ofstAux+1];
				//*t_EzRes = arEZ[ofstAux]; *(t_EzRes+1) = arEZ[ofstAux+1];
				//end test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				t_ExRes += 2;
				t_EzRes += 2;

				phEn += pWfr->eStep;
			}
			x += pWfr->xStep;
		}
		z += pWfr->zStep;
**/
	}

	return 0;
}
