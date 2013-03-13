/************************************************************************//**
 * File: srwlclient.cpp
 * Description: Demo C/C++ client 
 * Project: Synchrotron Radiation Workshop Library (SRWLib)
 * First release: October 2010
 *
 * SRW is Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * SRW C/C++ API (SRWLIB) is Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, G.Geloni, L.Samoylova
 * @version 0.04
 ***************************************************************************/

#include "srwlib.h"
#include <iostream>
#include <fstream>
#include <cstring> //necessary for strcpy, etc...
#include <cstdlib> //necessary for atoi
#include <sstream>
#include <iomanip>
using namespace std;

/************************************************************************//**
 * Auxiliary function dedicated to process errors reported by Library
 ***************************************************************************/
void ProcRes(int er)
{
	char ErrorBuf[2048];

	if(er == 0) return;
	else
	{
		srwlUtiGetErrText(ErrorBuf, er);

		cout << endl;
		if(er < 0) 
		{//Print Warning:
			cout << "WARNING: " << ErrorBuf << endl;
		}
		else 
		{//Just print Error Message:
			cout << "ERROR: " << ErrorBuf << endl;
			//cout << "Press Enter to exit" << endl;
			//getchar();
			//exit(0);
		}
	}
}

/************************************************************************//**
 * Wavefront modification (re-allocation) function; to be called by pointer from SRWLIB
 ***************************************************************************/
int ModifySRWLWfr(int action, SRWLWfr* pWfr, char pol)
{
	if(pWfr == 0) return -1; //returning non-zero means Wfr modification did not succeed; no throwing allowed here
	if((action < 0) || (action > 2)) return -1;

	long numTot = pWfr->mesh.ne*pWfr->mesh.nx*pWfr->mesh.ny*2;
	if(numTot <= 0) return 0; //or delete the previous array, still?

	int ExNeeded = ((pol == 0) || (pol == 'x') || (pol == 'X'))? 1 : 0;
	int EyNeeded = ((pol == 0) || (pol == 'y') || (pol == 'Y') || (pol == 'z') || (pol == 'Z'))? 1 : 0;

		//DEBUG
		cout << "nx= " << pWfr->mesh.nx << ", ny= " << pWfr->mesh.ny << endl; 

	if(action == 0) 
	{//just delete existing wavefront data
		if(ExNeeded)
		{
			if(pWfr->arEx) delete[] pWfr->arEx;
			pWfr->arEx = 0;
			if(pWfr->arMomX) delete[] pWfr->arMomX;
			pWfr->arMomX = 0;
		}
		if(EyNeeded)
		{
			if(pWfr->arEy) delete[] pWfr->arEy;
			pWfr->arEy = 0;
			if(pWfr->arMomY) delete[] pWfr->arMomY;
			pWfr->arMomY = 0;
		}
	}
	else if(action == 1)
	{//allocate new wavefront data (without checking/deleting any existing data)
		if(ExNeeded)
		{
			pWfr->arEx = (char*)(new float[numTot]);
			pWfr->arMomX = new double[11*pWfr->mesh.ne];
		}
		if(EyNeeded)
		{
			pWfr->arEy = (char*)(new float[numTot]);
			pWfr->arMomY = new double[11*pWfr->mesh.ne];
		}
	}
	else if(action == 2)
	{//modify wavefront size (numbers of points vs photon energy, horizontal or vertical position)
		if(ExNeeded)
		{//using realloc could perhaps be more efficient here
			if(pWfr->arEx) delete[] pWfr->arEx;
			pWfr->arEx = (char*)(new float[numTot]);
			if(pWfr->arMomX) delete[] pWfr->arMomX;
			pWfr->arMomX = new double[11*pWfr->mesh.ne];
		}
		if(EyNeeded)
		{//using realloc could perhaps be more efficient here
			if(pWfr->arEy) delete[] pWfr->arEy;
			pWfr->arEy = (char*)(new float[numTot]);
			if(pWfr->arMomY) delete[] pWfr->arMomY;
			pWfr->arMomY = new double[11*pWfr->mesh.ne];
		}
	}
	return 0;
}

/************************************************************************//**
 * Auxiliary function to read line from file and to extract a number from it
 ***************************************************************************/
template<class T> void AuxReadLineAndExtractNumber(T& numOut, ifstream& fIn, const string& sCom) //throw(...) 
{
	string sIn;
	getline(fIn, sIn); if(!fIn.good()) throw -1;

	size_t posComStart = sIn.find(sCom);
	if(posComStart != 0) throw -1;

	posComStart = sCom.length();
	size_t posComEnd = sIn.find(sCom, posComStart);

	istringstream is(sIn.substr(posComStart, posComEnd - 1));
	is >> numOut;
}

/************************************************************************//**
 * Auxiliary function to read 3D magnetic field data from ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxReadInMagFld3D(SRWLMagFld3D* pFld, const char* strFileName)
{
	if((strFileName == 0) || (pFld == 0)) return -1;
		//cout << strFileName << endl;
	ifstream f(strFileName);
	if(!f.is_open()) return -1;

	string sRead;
	double xStart = 0, xStep = 0, yStart = 0, yStep = 0, zStart = 0, zStep = 0;
	int xNp = 1, yNp = 1, zNp = 1;
	pFld->arBx = 0; pFld->arBy = 0; pFld->arBx = 0;

	try
	{
		getline(f, sRead); if(!f.good()) return -1; //1st line: just pass 
		AuxReadLineAndExtractNumber(xStart, f, "#"); //2nd line: initial X position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(xStep, f, "#"); //3rd line: step vs X [m]
		AuxReadLineAndExtractNumber(xNp, f, "#"); //4th line: number of points vs X

		AuxReadLineAndExtractNumber(yStart, f, "#"); //5th line: initial Y position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(yStep, f, "#"); //6th line: step vs Y [m]
		AuxReadLineAndExtractNumber(yNp, f, "#"); //7th line: number of points vs Y

		AuxReadLineAndExtractNumber(zStart, f, "#"); //8th line: initial Z position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(zStep, f, "#"); //9th line: step vs Z [m]
		AuxReadLineAndExtractNumber(zNp, f, "#"); //10th line: number of points vs Z

		pFld->rx = xStep*(xNp - 1);
		pFld->nx = xNp;
		pFld->ry = yStep*(yNp - 1);
		pFld->ny = yNp;
		pFld->rz = zStep*(zNp - 1);
		pFld->nz = zNp;

		long bNp = xNp*yNp*zNp;
		pFld->arBx = new double[bNp];
		pFld->arBy = new double[bNp];
		pFld->arBz = new double[bNp];

		double *t_arBx = pFld->arBx, *t_arBy = pFld->arBy, *t_arBz = pFld->arBz;
		for(long i=0; i<bNp; i++)
		{
			//lines from 11th: Magnetic Field Components Bx, By, Bz
			getline(f, sRead); if(!f.good()) return -1;
			istringstream is(sRead);
			is >> *(t_arBx++); 
			is >> *(t_arBy++); 
			is >> *(t_arBz++); 
		}
		f.close();
	}
	catch(int erNo) 
	{
		if(pFld->arBx != 0) { delete[] pFld->arBx; pFld->arBx = 0;}
		if(pFld->arBy != 0) { delete[] pFld->arBy; pFld->arBy = 0;}
		if(pFld->arBz != 0) { delete[] pFld->arBz; pFld->arBz = 0;}
		if(f.is_open()) f.close();
		return erNo;
	}
	return 0;
}

/************************************************************************//**
 * Auxiliary function to save trajectory data to ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxSaveTrajData(SRWLPrtTrj* pTrj, const char* strFileName)
{
	if((strFileName == 0) || (pTrj == 0)) return -1;

	ofstream f(strFileName);
	if(!f.is_open()) return -1;

	int nCh = 12; //number of characters in each output value
	f.precision(nCh);

	f << "#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [m]" << endl;

	double ctStep = (pTrj->np > 0)? (pTrj->ctEnd - pTrj->ctStart)/(pTrj->np - 1) : 0;
	double ct = pTrj->ctStart;
	double *t_arX = pTrj->arX, *t_arXp = pTrj->arXp;
	double *t_arY = pTrj->arY, *t_arYp = pTrj->arYp;
	double *t_arZ = pTrj->arZ, *t_arZp = pTrj->arZp;
	for(long i=0; i<pTrj->np; i++)
	{
		f << ct << '\t' << *(t_arX++) << '\t' << *(t_arXp++) << '\t' << *(t_arY++) << '\t' << *(t_arYp++) << '\t' << *(t_arZ++) << '\t' << *(t_arZp++) << endl;
		ct += ctStep;
	}
	f << ends;

	if(f.is_open()) f.close();
	return 0;
}

/************************************************************************//**
 * Auxiliary function to save intensity data to ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxSaveIntensData(float* arI, double eSt, double eFi, int ne, double xSt, double xFi, int nx, double ySt, double yFi, int ny, const char* strFileName)
{
	if((strFileName == 0) || (arI == 0)) return -1;

	ofstream f(strFileName);
	if(!f.is_open()) return -1;

	int nCh = 8; //number of characters in each output value
	f.precision(nCh);
	f << "C-aligned Intensity (inner loop is vs photon energy, outer loop vs vertical position)" << endl;
	f << '#' << eSt << " #Initial Photon Energy [eV]\n";
    f << '#' << eFi << " #Final Photon Energy [eV]\n";
    f << '#' << ne << " #Number of points vs Photon Energy\n";
    f << '#' << xSt << " #Initial Horizontal Position [m]\n";
    f << '#' << xFi << " #Final Horizontal Position [m]\n";
    f << '#' << nx << " #Number of points vs Horizontal Position\n";
    f << '#' << ySt << " #Initial Vertical Position [m]\n";
    f << '#' << yFi << " #Final Vertical Position [m]\n";
    f << '#' << ny << " #Number of points vs Vertical Position\n";

	float *t_atI = arI;
	for(long i=0; i<(ne*nx*ny); i++) f << " " << *(t_atI++) << endl;
	f << ends;
	if(f.is_open()) f.close();
	return 0;
}

/************************************************************************//**
 * Example#1: Calculating electron trajectory in 3D magnetic field of an APPLE-II undulator
 ***************************************************************************/
int SRWLIB_Example01(const char* strFolder)
{
	cout << "SRWLIB C Client Example # 1:" << endl; 
	cout << "Calculating electron trajectory in 3D magnetic field of an APPLE-II undulator" << endl; 

	//**********************Input Parameters:
	char strExampleFolderName[] = "data_example_01"; //example sub-folder name
	//char *arFldInFileNames[] = {"epu49term1.dat", "epu49cen.dat", "epu49term2.dat"}; //3D Magnetic Field data file names
	char *arFldInFileNames[] = {"epu49HEtot.dat"}; //3D Magnetic Field data file names
	char strTrajOutFileName[] = "ex01_res_traj.dat"; //file name for output trajectory data

	int numPer = 40; //Number of ID Periods
	double xcID = 0.; //[m] Transverse Coordinates of ID Center
	double ycID = 0.;
	double zcID = 0.; //[m] Longitudinal Coordinate of ID Center

	SRWLPrtTrj partTraj; //Trajectory structure (to be used both for input and output)
	partTraj.np = 10001; //Number of Points for Trajectory calculation
	
	int fieldInterpMeth = 3; //Magnetic Field Interpolation Method, to be entered into 3D field structures below (to be used e.g. for trajectory calculation):
	//1- bi-linear (3D), 2- (bi-)quadratic (3D), 3- (bi-)cubic (3D)
	double arPrecPar[2]; //General Precision parameters for Trajectory calculation:
	arPrecPar[0] = 1; //Number of precision parameters
	arPrecPar[1] = 1; //Integration method No:
		//1- fourth-order Runge-Kutta (precision is driven by number of points)
		//2- fifth-order Runge-Kutta
	//arPrecPar[1],[2],[3],[4],[5]: absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] (yet to be tested!!) - to be taken into account only for R-K fifth order or higher
	//arPrecPar[6]: tolerance (default = 1) for R-K fifth order or higher
	//arPrecPar[7]: max. number of auto-steps for R-K fifth order or higher (default = 5000)
	
	SRWLParticle &part = partTraj.partInitCond; //Particle structure for Initial Conditions
	//Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on):
	part.x = 0.004; //[m]
	part.y = 0.001; //[m]
	//Initial Transverse Velocities:
	part.xp = 0.;
	part.yp = 0.;
	part.gamma = 3./0.51099890221e-03; //Relative Energy
	part.relE0 = 1; //Electron rest mass
	part.nq = -1; //Electron charge
	partTraj.ctStart = 0.; //Start Time
	//partTraj.ctEnd will be defined later on, after ID length will be calculated (particle is assumed to pass through all ID)
	//**********************End of Input Parameters

	char strBuf[2048];
	char strSep[] = "/\0";
	strcpy(strBuf, strFolder);
	strcat(strBuf, strExampleFolderName);
	size_t lenStrFolder = strlen(strFolder);
	const char *pSep = (lenStrFolder > 0)? strFolder + (lenStrFolder - 1) : strSep;
	strcat(strBuf, pSep);
	char *pFileName = strBuf + strlen(strBuf);

	//Creating array of 3 pieces of 3D Magnetic Field
	//SRWLMagFld3D arMagFld[3];
	//for(int i=0; i<3; i++)
	//{//Setting all pointers to 0 before any allocation - to allow for selective deallocation in case of error
	//	SRWLMagFld3D *pMagFld = arMagFld + i;
	//	pMagFld->arBx = pMagFld->arBy = pMagFld->arBz = 0;
	//	pMagFld->arX = pMagFld->arY = pMagFld->arZ = 0; //if these pointers are not zero, regular mesh data will be ignored
	//	pMagFld->interp = fieldInterpMeth;
	//}

	//Allocating Trajectory arrays
	long totNumTrajPoints = partTraj.np*6;
	double *arTrajData = new double[totNumTrajPoints]; //one big array instead of 6 smaller ones
	double *t_arTrajData = arTrajData;
	//Setting pointers
	partTraj.arX = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arXp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arY = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arYp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZ = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZp = t_arTrajData;

	char arFldTypes[3];
	double arXcID[3], arYcID[3], arZcID[3];
	void *arPtrVoid[3]; //auxiliary array of pointers for defining Container
	//SRWLMagFld3D *t_arMagFld = arMagFld;

	SRWLMagFld3D arMagFld[1];

	cout << "   Reading magnetic field data from files ... "; 
	//for(int j=0; j<3; j++)
	for(int j=0; j<1; j++)
	{
		SRWLMagFld3D *pMagFld = arMagFld + j;
		pMagFld->arBx = pMagFld->arBy = pMagFld->arBz = 0; //Setting all pointers to 0 before any allocation - to allow for selective deallocation in case of error
		pMagFld->arX = pMagFld->arY = pMagFld->arZ = 0; //if these pointers are not zero, regular mesh data will be ignored
		pMagFld->interp = fieldInterpMeth;

		//Preparing full path string (strBuf)
		strcpy(pFileName, arFldInFileNames[j]);
		//Reading field data from file:
		if(AuxReadInMagFld3D(pMagFld, strBuf)) 
			cout << "Error reading 3D magnetic field data from file" << endl;

		arPtrVoid[j] = (void*)pMagFld;
		arFldTypes[j] = 'a'; //all field types are "Arbitrary 3D"
		arXcID[j] = xcID;
		arYcID[j] = ycID;
		arZcID[j] = zcID;
	}
	cout << "done" << endl; 

	arMagFld[0].nRep = 1; //1st termination
	//arMagFld[1].nRep = numPer; //central part, repeated numPer times
	//arMagFld[2].nRep = 1; //2nd termination
	//double per = arMagFld[1].rz;

	//Defining Magnetic Field Container
	SRWLMagFldC magFldCnt;
	magFldCnt.arMagFld = arPtrVoid; 
	magFldCnt.arMagFldTypes = arFldTypes;
	magFldCnt.arXc = arXcID;
	magFldCnt.arYc = arYcID;
	magFldCnt.arZc = arZcID;
	//double arMagZc[] = {zcID - 0.5*numPer*per - 0.5*(arMagFld[0].rz), zcID, zcID + 0.5*numPer*per + 0.5*(arMagFld[2].rz)};
	//magFldCnt.arZc = arMagZc;
	magFldCnt.nElem = 1; //3;

	//Defining remaining Initial Conditions based on ID length
	//part.z = zcID - (0.5*numPer + 1)*per - arMagFld[0].rz; //Initial Longitudinal Coordinate
	//partTraj.ctEnd = (numPer + 2)*per + arMagFld[0].rz + arMagFld[2].rz; //End Time
	part.z = zcID - 0.5*arMagFld[0].rz; //Initial Longitudinal Coordinate
	partTraj.ctEnd = arMagFld[0].rz;

	//SRWLIB function call (to run calculations)
	cout << "   Performing calculation ... "; 
	ProcRes(srwlCalcPartTraj(&partTraj, &magFldCnt, arPrecPar));
	cout << "done" << endl; 

	//Preparing full path to output ASCII file and Saving Trajectory data into it:
	strcpy(pFileName, strTrajOutFileName);
	cout << "   Saving trajectory data to a file ... "; 
	if(AuxSaveTrajData(&partTraj, strBuf)) 
		cout << "Error saving resulting trajectory data to file" << endl;
	cout << "done" << endl; 

	//**********************Deallocating memory
	//for(int k=0; k<3; k++)
	for(int k=0; k<1; k++)
	{
		SRWLMagFld3D *pMagFld = arMagFld + k;
		if(pMagFld->arBx != 0) delete[] pMagFld->arBx;
		if(pMagFld->arBy != 0) delete[] pMagFld->arBy;
		if(pMagFld->arBz != 0) delete[] pMagFld->arBz;
	}
	if(arTrajData != 0) delete[] arTrajData;
	return 0;
}

/************************************************************************//**
 * Example#2: Calculating electron trajectory in magnetic field of a segmented planar undulator with FODO lattice
 ***************************************************************************/
int SRWLIB_Example02(const char* strFolder)
{
	cout << "SRWLIB C Client Example # 2:" << endl; 
	cout << "Calculating electron trajectory in magnetic field of a segmented planar undulator with FODO lattice" << endl; 

	//**********************Input Parameters:
	char strExampleFolderName[] = "data_example_02"; //example data sub-folder name
	char strTrajOutFileName[] = "ex02_res_traj.dat"; //file name for output trajectory data
	
	int numSegm = 5; //Number of ID Segments
	int numPer = 100; //Number of Periods in one Segment (without counting for terminations)
	double undPer = 0.02; //Period Length [m]
	double xcID = 0; //Transverse Coordinates of ID Center [m]
	double ycID = 0;
	double zcID = 0; //Longitudinal Coordinate of ID Center [m]
	
	SRWLParticle part;
	part.x = 0.0001; //Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
	part.y = 0.0001;
	part.xp = 0; //Initial Transverse Velocities
	part.yp = 0;
	part.gamma = 3/0.51099890221e-03; //Relative Energy
	part.relE0 = 1; //Electron Rest Mass
	part.nq = -1; //Electron Charge

	int npTraj = 20001; //Number of Points for Trajectory calculation

	double arPrecPar[2]; //General Precision parameters for Trajectory calculation:
	arPrecPar[0] = 1; //Number of precision parameters
	arPrecPar[1] = 1; //Integration method No:
		//1- fourth-order Runge-Kutta (precision is driven by number of points)
		//2- fifth-order Runge-Kutta
	//arPrecPar[1],[2],[3],[4],[5]: absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] (yet to be tested!!) - to be taken into account only for R-K fifth order or higher
	//arPrecPar[6]: tolerance (default = 1) for R-K fifth order or higher
	//arPrecPar[7]: max. number of auto-steps for R-K fifth order or higher (default = 5000)

	SRWLMagFldH harm;
	harm.n = 1; //harmonic number
	harm.h_or_v = 'v'; //magnetic field plane: horzontal ('h') or vertical ('v')
	harm.B = 1.05; //magnetic field amplitude [T]
	harm.ph = 0; //phase [rad]
	harm.s = 1; //symmetry vs longitudinal position: 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
	harm.a = 1; //coefficient for transverse depenednce: B*cosh(2*Pi*n*a*y/per)*cos(2*Pi*n*z/per + ph)

	SRWLMagFldU und;
	und.arHarm = &harm; //array of field harmonics
	und.nHarm = 1; //number of field harmonics
	und.per = undPer; //period length [m]
	und.nPer = numPer; //number of periods

	SRWLMagFldM qf, qd;
	qf.G = 0.5; //field parameter: [T] for dipole, [T/m] for quadrupole (negative means defocusing for x), [T/m^2] for sextupole, [T/m^3] for octupole
	qd.G = -0.5; //field parameter: [T] for dipole, [T/m] for quadrupole (negative means defocusing for x), [T/m^2] for sextupole, [T/m^3] for octupole
	qf.m = qd.m = 2; //multipole order: 1 for dipole, 2 for quadrupoole, 3 for sextupole, 4 for octupole
	qf.n_or_s = qd.n_or_s = 'n'; //normal ('n') or skew ('s')
	qf.Leff = qd.Leff = 0.2; //effective length [m]
	qf.Ledge = qd.Ledge = 0; //"soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed

	double arZero[] = {0,0,0,0,0,0,0,0,0};
	double undLen = (numPer + 2)*undPer;
	double distBwSegm = 0.4; //Distance between Undulator Segments
	double undLenExt = undLen + distBwSegm;
	double arZcen[] = {-2*undLenExt, -1.5*undLenExt, -undLenExt, -0.5*undLenExt, 0, 0.5*undLenExt, undLenExt, 1.5*undLenExt, 2*undLenExt};

	SRWLMagFldC magFldCnt;
	void *vArMagFld[] = {(void*)(&und), (void*)(&qf), (void*)(&und), (void*)(&qd), (void*)(&und), (void*)(&qf), (void*)(&und), (void*)(&qd), (void*)(&und)};
	magFldCnt.arMagFld = vArMagFld; //array of pointers to magnetic field elements
	magFldCnt.arMagFldTypes = "umumumumu"; //types of magnetic field elements in arMagFld array
	magFldCnt.arXc = arZero; //horizontal center positions of magnetic field elements in arMagFld array
	magFldCnt.arYc = arZero; //vertical center positions of magnetic field elements in arMagFld array
	magFldCnt.arZc = arZcen; //longitudinal center positions of magnetic field elements in arMagFld array
	magFldCnt.nElem = 9; //number of magnetic field elements in arMagFld array

	part.z = arZcen[0] - 0.5*undLenExt; //Initial Longitudinal Coordinate (set before the ID)

	SRWLPrtTrj partTraj;
	partTraj.partInitCond = part; //particle type and initial conditions for which the trajectory should be (/is) calculated
	partTraj.np = npTraj; //number of trajectory points
	partTraj.ctStart = 0; //"Start Time" (c*t) for the calculation (0 corresponds to the time moment for which the initial conditions are defined)
	partTraj.ctEnd = partTraj.ctStart + 5*undLenExt; //End Time
	//Allocating Trajectory arrays
	long totNumTrajPoints = npTraj*6;
	double *arTrajData = new double[totNumTrajPoints]; //one big array instead of 6 smaller ones
	double *t_arTrajData = arTrajData;
	//Setting pointers //arrays of horizontal, vertical and longitudinal positions and relative velocities
	partTraj.arX = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arXp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arY = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arYp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZ = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZp = t_arTrajData;

	//SRWLIB function call (to run calculations)
	cout << "   Performing calculation ... "; 
	ProcRes(srwlCalcPartTraj(&partTraj, &magFldCnt, arPrecPar));
	cout << "done" << endl; 

	//Preparing full path to output ASCII file and Saving Trajectory data into it:
	char strBuf[2048];
	char strSep[] = "/\0";
	strcpy(strBuf, strFolder);
	strcat(strBuf, strExampleFolderName);
	size_t lenStrFolder = strlen(strFolder);
	const char *pSep = (lenStrFolder > 0)? strFolder + (lenStrFolder - 1) : strSep;
	strcat(strBuf, pSep);

	char *pFileName = strBuf + strlen(strBuf);
	strcpy(pFileName, strTrajOutFileName);
	cout << "   Saving trajectory data to a file ... "; 
	if(AuxSaveTrajData(&partTraj, strBuf)) cout << "Error saving resulting trajectory data to file" << endl;
	cout << "done" << endl; 

	//**********************Deallocating memory
	if(arTrajData != 0) delete[] arTrajData;
	return 0;
}

/************************************************************************//**
 * Example#3: Calculating synchrotron (undulator) radiation emitted by an electron travelling in ellipsoidal undulator
 ***************************************************************************/
int SRWLIB_Example03(const char* strFolder)
{
	cout << "SRWLIB C Client Example # 3:" << endl; 
	cout << "Calculating synchrotron (undulator) radiation emitted by an electron travelling in ellipsoidal undulator" << endl; 

	//**********************Input Parameters:
	char strExampleFolderName[] = "data_example_03"; //example data sub-folder name
	char strIntOutFileName1[] = "ex03_res_int1.dat"; //file name for output SR intensity data
	char strIntOutFileName2[] = "ex03_res_int2.dat"; //file name for output SR intensity data
	
	//***********Undulator
	int numPer = 40; //Number of ID Periods (without counting for terminations
	double undPer = 0.049; //Period Length [m]
	double Bx = 0.57/3.; //Peak Horizontal field [T]
	double By = 0.57; //Peak Vertical field [T]
	double phBx = 0; //Initial Phase of the Horizontal field component
	double phBy = 0; //Initial Phase of the Vertical field component
	int sBx = -1; //Symmetry of the Horizontal field component vs Longitudinal position
	int sBy = 1; //Symmetry of the Vertical field component vs Longitudinal position
	double xcID = 0; //Transverse Coordinates of Undulator Center [m]
	double ycID = 0;
	double zcID = 0; //Longitudinal Coordinate of Undulator Center [m]

	SRWLMagFldH harmY, harmX;
	harmY.n = harmX.n = 1; //harmonic number
	harmY.h_or_v = 'v'; harmX.h_or_v = 'h'; //magnetic field plane: horzontal ('h') or vertical ('v')
	harmY.B = By; harmX.B = Bx; //magnetic field amplitude [T]
	harmY.ph = harmX.ph = 0; //phase [rad]
	harmY.s = sBy; harmX.s = sBx; //symmetry vs longitudinal position: 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
	harmY.a = harmX.a = 1; //coefficient for transverse depenednce: B*cosh(2*Pi*n*a*y/per)*cos(2*Pi*n*z/per + ph)
	SRWLMagFldH arH[] = {harmY, harmX};

	SRWLMagFldU und; //Ellipsoidal Undulator
	und.arHarm = arH; //arH; //arHarmonics; //array of field harmonics
	und.nHarm = 2; //number of field harmonics
	und.per = undPer; //period length [m]
	und.nPer = numPer; //number of periods

	SRWLMagFldC magFldCnt;
	void *vArMagFld[] = {(void*)(&und)};
	magFldCnt.arMagFld = vArMagFld; //array of pointers to magnetic field elements
	magFldCnt.arMagFldTypes = "u"; //types of magnetic field elements in arMagFld array
	double auxArXcID[] = {xcID};
	double auxArYcID[] = {ycID};
	double auxArZcID[] = {zcID};
	magFldCnt.arXc = auxArXcID; //horizontal center positions of magnetic field elements in arMagFld array
	magFldCnt.arYc = auxArYcID; //vertical center positions of magnetic field elements in arMagFld array
	magFldCnt.arZc = auxArZcID; //longitudinal center positions of magnetic field elements in arMagFld array
	magFldCnt.nElem = 1; //number of magnetic field elements in arMagFld array

	//***********Electron Beam
	SRWLPartBeam elecBeam;
	elecBeam.Iavg = 0.5; //Average Current [A]
	elecBeam.partStatMom1.x = 0.; //Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
	elecBeam.partStatMom1.y = 0.;
	elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 4); //Initial Longitudinal Coordinate (set before the ID)
	elecBeam.partStatMom1.xp = 0; //Initial Relative Transverse Velocities
	elecBeam.partStatMom1.yp = 0;
	elecBeam.partStatMom1.gamma = 3./0.51099890221e-03; //Relative Energy
	elecBeam.partStatMom1.relE0 = 1; //Rest mass (energy) in units of electron rest mass: =1 for electron, =1836.1526988 (=938.272013/0.510998902) for proton
	elecBeam.partStatMom1.nq = -1; //Charge of the particle related to absolute value of electron charge: =-1 for electron, =1 for positron and for proton

	//***********Precision
	int meth = 1; //SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
	double relPrec = 0.01; //relative precision
	double zStartInteg = 0; //longitudinal position to start integration (effective if < zEndInteg)
	double zEndInteg = 0; //longitudinal position to finish integration (effective if > zStartInteg)
	int npTraj = 20000;
	double sampFactNxNyForProp = 0; //sampling factor for adjusting nx, ny (effective if > 0)
	double arPrecPar[] = {meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp};

	//***********Wavefronts
	SRWLWfr wfr1; //For spectrum vs photon energy
	wfr1.mesh.ne = 10000; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr1.mesh.nx = wfr1.mesh.ny = 1;
	wfr1.mesh.zStart = 20.; //Longitudinal Position [m] at which SR has to be calculated
	wfr1.mesh.eStart = 10.; //Initial Photon Energy [eV]
	wfr1.mesh.eFin = 3000.; //Final Photon Energy [eV]
	wfr1.mesh.xStart = 0.; //Initial Horizontal Position [m]
	wfr1.mesh.xFin = 0; //Final Horizontal Position [m]
	wfr1.mesh.yStart = 0; //Initial Vertical Position [m]
	wfr1.mesh.yFin = 0; //Final Vertical Position [m]
	wfr1.partBeam = elecBeam;
	wfr1.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr1.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	long numTot = wfr1.mesh.ne*wfr1.mesh.nx*wfr1.mesh.ny*2;
	float *arEx1 = new float[numTot];
	float *arEy1 = new float[numTot];
	wfr1.arEx = (char*)(arEx1); //horizontal and vertical electric field component arrays
	wfr1.arEy = (char*)(arEy1);
	wfr1.arElecPropMatr = new double[20];
	wfr1.arWfrAuxData = new double[30];
	wfr1.arMomX = new double[11*wfr1.mesh.ne];
	wfr1.arMomY = new double[11*wfr1.mesh.ne];

	SRWLWfr wfr2; //For intensity distribution at fixed photon energy
	wfr2.mesh.ne = 1; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr2.mesh.nx = wfr2.mesh.ny = 101;
	wfr2.mesh.zStart = 20.; //Longitudinal Position [m] at which SR has to be calculated
	wfr2.mesh.eStart = 1090.; //Initial Photon Energy [eV]
	wfr2.mesh.eFin = 1090.; //Final Photon Energy [eV]
	wfr2.mesh.xStart = -0.001; //Initial Horizontal Position [m]
	wfr2.mesh.xFin = 0.001; //Final Horizontal Position [m]
	wfr2.mesh.yStart = -0.001; //Initial Vertical Position [m]
	wfr2.mesh.yFin = 0.001; //Final Vertical Position [m]
	wfr2.partBeam = elecBeam;
	wfr2.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr2.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	numTot = wfr2.mesh.ne*wfr2.mesh.nx*wfr2.mesh.ny*2;
	float *arEx2 = new float[numTot];
	float *arEy2 = new float[numTot];
	wfr2.arEx = (char*)(arEx2); //horizontal and vertical electric field component arrays
	wfr2.arEy = (char*)(arEy2);
	wfr2.arElecPropMatr = new double[20];
	wfr2.arWfrAuxData = new double[30];
	wfr2.arMomX = new double[11*wfr2.mesh.ne];
	wfr2.arMomY = new double[11*wfr2.mesh.ne];

	//**********************Calculation (SRWLIB function calls)
	cout << "   Performing Electric Field calculation ... ";
	ProcRes(srwlCalcElecFieldSR(&wfr1, 0, &magFldCnt, arPrecPar));
	cout << "done" << endl;
	cout << "   Extracting Intensity from calculated Electric Field ... ";
	float *arI1 = new float[wfr1.mesh.ne];
	ProcRes(srwlCalcIntFromElecField((char*)arI1, &wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart));
	cout << "done" << endl;

	cout << "   Performing Electric Field calculation ... ";
	//arPrecPar[6] = 1.3; //sampling factor for adjusting nx, ny (effective if > 0)
	ProcRes(srwlCalcElecFieldSR(&wfr2, 0, &magFldCnt, arPrecPar));
	cout << "done" << endl;
	cout << "   Extracting Intensity from calculated Electric Field ... ";
	float *arI2 = new float[wfr2.mesh.nx*wfr2.mesh.ny]; //"flat" array to take 2D intensity data
	ProcRes(srwlCalcIntFromElecField((char*)arI2, &wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0));
	cout << "done" << endl;

	//**********************Preparing full path to output ASCII file and Saving Intensity data into it:
	char strBuf[2048];
	char strSep[] = "/\0";
	strcpy(strBuf, strFolder);
	strcat(strBuf, strExampleFolderName);
	size_t lenStrFolder = strlen(strFolder);
	const char *pSep = (lenStrFolder > 0)? strFolder + (lenStrFolder - 1) : strSep;
	strcat(strBuf, pSep);

	char *pFileName = strBuf + strlen(strBuf);
	strcpy(pFileName, strIntOutFileName1);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI1, wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne, wfr1.mesh.xStart, wfr1.mesh.xFin, wfr1.mesh.nx, wfr1.mesh.yStart, wfr1.mesh.yFin, wfr1.mesh.ny, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl;

	strcpy(pFileName, strIntOutFileName2);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI2, wfr2.mesh.eStart, wfr2.mesh.eFin, wfr2.mesh.ne, wfr2.mesh.xStart, wfr2.mesh.xFin, wfr2.mesh.nx, wfr2.mesh.yStart, wfr2.mesh.yFin, wfr2.mesh.ny, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl; 

	//**********************Deallocating memory
	delete[] arEx1;
	delete[] arEy1;
	delete[] wfr1.arElecPropMatr;
	delete[] wfr1.arWfrAuxData;
	delete[] wfr1.arMomX;
	delete[] wfr1.arMomY;
	delete[] arI1;
	delete[] arEx2;
	delete[] arEy2;
	delete[] wfr2.arElecPropMatr;
	delete[] wfr2.arWfrAuxData;
	delete[] wfr2.arMomX;
	delete[] wfr2.arMomY;
	delete[] arI2;
	return 0;
}

/************************************************************************//**
 * Example#4: Calculating synchrotron (undulator) radiation electric field (from one electron)
 * and simulating wavefront propagation through a simple optical system
 ***************************************************************************/
int SRWLIB_Example04(const char* strFolder)
{
	cout << "SRWLIB C Client Example # 4:" << endl; 
	cout << "Calculating synchrotron (undulator) radiation electric field (from one electron) and performing simulation of wavefront propagation through a simple optical system" << endl; 

	//**********************Input Parameters:
	char strExampleFolderName[] = "data_example_04"; //example data sub-folder name
	char strIntOutFileName1[] = "ex04_res_int1.dat"; //file name for output SR intensity data
	char strIntOutFileName2[] = "ex04_res_int2.dat"; //file name for output SR intensity data

	//***********Undulator
	int numPer = 40; //Number of ID Periods (without counting for terminations
	double undPer = 0.049; //Period Length [m]
	double Bx = 0.57/3.; //Peak Horizontal field [T]
	double By = 0.57; //Peak Vertical field [T]
	double phBx = 0; //Initial Phase of the Horizontal field component
	double phBy = 0; //Initial Phase of the Vertical field component
	int sBx = -1; //Symmetry of the Horizontal field component vs Longitudinal position
	int sBy = 1; //Symmetry of the Vertical field component vs Longitudinal position
	double xcID = 0; //Transverse Coordinates of Undulator Center [m]
	double ycID = 0;
	double zcID = 0; //Longitudinal Coordinate of Undulator Center [m]

	SRWLMagFldH harmY, harmX;
	harmY.n = harmX.n = 1; //harmonic number
	harmY.h_or_v = 'v'; harmX.h_or_v = 'h'; //magnetic field plane: horzontal ('h') or vertical ('v')
	harmY.B = By; harmX.B = Bx; //magnetic field amplitude [T]
	harmY.ph = harmX.ph = 0; //phase [rad]
	harmY.s = sBy; harmX.s = sBx; //symmetry vs longitudinal position: 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
	harmY.a = harmX.a = 1; //coefficient for transverse depenednce: B*cosh(2*Pi*n*a*y/per)*cos(2*Pi*n*z/per + ph)
	SRWLMagFldH arH[] = {harmY, harmX};

	SRWLMagFldU und; //Ellipsoidal Undulator
	und.arHarm = arH; //arH; //arHarmonics; //array of field harmonics
	und.nHarm = 2; //number of field harmonics
	und.per = undPer; //period length [m]
	und.nPer = numPer; //number of periods

	SRWLMagFldC magFldCnt;
	void *vArMagFld[] = {(void*)(&und)};
	magFldCnt.arMagFld = vArMagFld; //array of pointers to magnetic field elements
	magFldCnt.arMagFldTypes = "u"; //types of magnetic field elements in arMagFld array
	double auxArXcID[] = {xcID};
	double auxArYcID[] = {ycID};
	double auxArZcID[] = {zcID};
	magFldCnt.arXc = auxArXcID; //horizontal center positions of magnetic field elements in arMagFld array
	magFldCnt.arYc = auxArYcID; //vertical center positions of magnetic field elements in arMagFld array
	magFldCnt.arZc = auxArZcID; //longitudinal center positions of magnetic field elements in arMagFld array
	magFldCnt.nElem = 1; //number of magnetic field elements in arMagFld array

	//***********Electron Beam
	SRWLPartBeam elecBeam;
	elecBeam.Iavg = 0.5; //Average Current [A]
	elecBeam.partStatMom1.x = 0.; //Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
	elecBeam.partStatMom1.y = 0.;
	elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 4); //Initial Longitudinal Coordinate (set before the ID)
	elecBeam.partStatMom1.xp = 0; //Initial Relative Transverse Velocities
	elecBeam.partStatMom1.yp = 0;
	elecBeam.partStatMom1.gamma = 3./0.51099890221e-03; //Relative Energy
	elecBeam.partStatMom1.relE0 = 1; //Rest mass (energy) in units of electron rest mass: =1 for electron, =1836.1526988 (=938.272013/0.510998902) for proton
	elecBeam.partStatMom1.nq = -1; //Charge of the particle related to absolute value of electron charge: =-1 for electron, =1 for positron and for proton

	//***********Precision for SR calculation
	int meth = 1; //SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
	double relPrec = 0.01; //relative precision
	double zStartInteg = 0; //longitudinal position to start integration (effective if < zEndInteg)
	double zEndInteg = 0; //longitudinal position to finish integration (effective if > zStartInteg)
	int npTraj = 20000;
	double sampFactNxNyForProp = 1; //sampling factor for adjusting nx, ny (effective if > 0)
	double arPrecPar[] = {meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp};

	//***********Initial Wavefront
	SRWLWfr wfr; //For spectrum vs photon energy
	wfr.mesh.ne = 1; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr.mesh.nx = wfr.mesh.ny = 100;
	wfr.mesh.zStart = 20.; //Longitudinal Position [m] at which SR has to be calculated
	wfr.mesh.eStart = 1090.; //Initial Photon Energy [eV]
	wfr.mesh.eFin = 1090.; //Final Photon Energy [eV]
	wfr.mesh.xStart = -0.001; //Initial Horizontal Position [m]
	wfr.mesh.xFin = 0.001; //Final Horizontal Position [m]
	wfr.mesh.yStart = -0.001; //Initial Vertical Position [m]
	wfr.mesh.yFin = 0.001; //Final Vertical Position [m]
	wfr.partBeam = elecBeam;
	wfr.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	long numTot = wfr.mesh.ne*wfr.mesh.nx*wfr.mesh.ny*2;
	//float *arEx1 = new float[numTot];
	//float *arEy1 = new float[numTot];
	//wfr.arEx = (char*)(arEx1); //horizontal and vertical electric field component arrays
	//wfr.arEy = (char*)(arEy1);
	wfr.arEx = (char*)(new float[numTot]); //horizontal and vertical electric field component arrays
	wfr.arEy = (char*)(new float[numTot]);
	wfr.arElecPropMatr = new double[20];
	wfr.arWfrAuxData = new double[30];
	wfr.arMomX = new double[11*wfr.mesh.ne];
	wfr.arMomY = new double[11*wfr.mesh.ne];

	//***********Optical Elements and Propagation Parameters
	SRWLOptL optLens; //Lens
	optLens.Fx = wfr.mesh.zStart/2.; //Lens focal lengths
	optLens.Fy = wfr.mesh.zStart/2.;
	optLens.x = 0; //Transverse coordinates of center
	optLens.y = 0;

	SRWLOptD optDrift; //Drift space
	optDrift.L = wfr.mesh.zStart; 

	double propagParLens[] = {1, 1, 1., 0, 0, 1., 1.5, 1., 1.5, 0, 0, 0};
	double propagParDrift[] = {1, 1, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0};
	//double propagParLens[] = {0, 0, 1., 1, 0, 1., 1.5, 1., 1.5, 0, 0, 0};
	//double propagParDrift[] = {0, 0, 1., 1, 0, 1., 1., 1., 1., 0, 0, 0};
	//Wavefront Propagation Parameters:
	//[0]: Auto-Resize (1) or not (0) Before propagation
	//[1]: Auto-Resize (1) or not (0) After propagation
	//[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
	//[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
	//[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
	//[5]: Horizontal Range modification factor at Resizing (1. means no modification)
	//[6]: Horizontal Resolution modification factor at Resizing
	//[7]: Vertical Range modification factor at Resizing
	//[8]: Vertical Resolution modification factor at Resizing
	//[9]: Type of wavefront Shift before Resizing (not yet implemented)
	//[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
	//[11]: New Vertical wavefront Center position after Shift (not yet implemented)
	double *arPropPar[] = {propagParLens, propagParDrift};

	SRWLOptC optBL; //Beamline (container)
	void *arOptEl[] = {(void*)(&optLens), (void*)(&optDrift)};
	optBL.arOpt = arOptEl; //array of pointers to optical elements
	char *arOptElT[] = {"lens", "drift"};
	optBL.arOptTypes = arOptElT; //array of types of optical elements (C strings) in arOpt array
	optBL.nElem = 2; //number of magnetic field elements in arMagFld array
	optBL.arProp = arPropPar; //array of arrays of propagation parameters to be used for individual optical elements
	optBL.nProp = 2; //number of propagation instructions (length of arProp array); 

	//**********************Calculation (SRWLIB function calls) and post-processing
	cout << "   Performing Initial Electric Field calculation ... ";
	ProcRes(srwlCalcElecFieldSR(&wfr, 0, &magFldCnt, arPrecPar));
	cout << "done" << endl;
	cout << "   Extracting Intensity from calculated Initial Electric Field ... ";
	int nx0 = wfr.mesh.nx, ny0 = wfr.mesh.ny;
	double xSt0 = wfr.mesh.xStart, xFi0 = wfr.mesh.xFin;
	double ySt0 = wfr.mesh.yStart, yFi0 = wfr.mesh.yFin;
	float *arI0 = new float[nx0*ny0]; //"flat" array to take 2D intensity data
	ProcRes(srwlCalcIntFromElecField((char*)arI0, &wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0));
	cout << "done" << endl;

	cout << "   Simulating Electric Field Wavefront Propagation ... ";
	ProcRes(srwlPropagElecField(&wfr, &optBL));
	cout << "done" << endl;
	cout << "   Extracting Intensity from propagated Electric Field and saving it to files ... ";
	float *arI1 = new float[wfr.mesh.nx*wfr.mesh.ny]; //"flat" array to take 2D intensity data
	ProcRes(srwlCalcIntFromElecField((char*)arI1, &wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0));
	cout << "done" << endl;

	//**********************Preparing full path to output ASCII file and Saving Intensity data into it:
	char strBuf[2048];
	char strSep[] = "/\0";
	strcpy(strBuf, strFolder);
	strcat(strBuf, strExampleFolderName);
	size_t lenStrFolder = strlen(strFolder);
	const char *pSep = (lenStrFolder > 0)? strFolder + (lenStrFolder - 1) : strSep;
	strcat(strBuf, pSep);

	char *pFileName = strBuf + strlen(strBuf);
	strcpy(pFileName, strIntOutFileName1);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI0, wfr.mesh.eStart, wfr.mesh.eFin, 1, xSt0, xFi0, nx0, ySt0, yFi0, ny0, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl;

	strcpy(pFileName, strIntOutFileName2);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI1, wfr.mesh.eStart, wfr.mesh.eFin, 1, wfr.mesh.xStart, wfr.mesh.xFin, wfr.mesh.nx, wfr.mesh.yStart, wfr.mesh.yFin, wfr.mesh.ny, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl; 

	//**********************Deallocating memory
	if(wfr.arEx != 0) delete[] wfr.arEx;
	if(wfr.arEy != 0) delete[] wfr.arEy;
	if(wfr.arWfrAuxData != 0) delete[] wfr.arWfrAuxData;
	if(wfr.arMomX != 0) delete[] wfr.arMomX;
	if(wfr.arMomY != 0) delete[] wfr.arMomY;
	delete[] arI0;
	delete[] arI1;
	return 0;
}

/************************************************************************//**
 * Example#5: Calculating electron trajectory and spontaneous emission
 * from a very long segmented undulator (transversely-uniform magnetic field defined)
 ***************************************************************************/
int SRWLIB_Example05(const char* strFolder)
{
	cout << "SRWLIB C Client Example # 5:" << endl; 
	cout << "Calculating electron trajectory and spontaneous emission from a very long segmented undulator (transversely-uniform magnetic field defined)" << endl; 

	//**********************Input Parameters:
	char strExampleFolderName[] = "data_example_05"; //example data sub-folder name
	char arFldInFileName[] = "segmented.dat"; //3D Magnetic Field data file names
	char strTrajOutFileName[] = "ex05_res_traj.dat"; //file name for output trajectory data
	char strIntOutFileName1[] = "ex05_res_int1.dat"; //file name for output SR intensity data
	char strIntOutFileName2[] = "ex05_res_int2.dat"; //file name for output SR intensity data

	//**********************Defining Magnetic Field:
	double xcID = 0; //Transverse Coordinates of ID Center [m]
	double ycID = 0;
	double zcID = 0; //Longitudinal Coordinate of ID Center [m]

	cout << "   Reading magnetic field data from file ... "; 
	//Preparing full path string (strBuf)
	char strBuf[2048];
	char strSep[] = "/\0";
	strcpy(strBuf, strFolder);
	strcat(strBuf, strExampleFolderName);
	size_t lenStrFolder = strlen(strFolder);
	const char *pSep = (lenStrFolder > 0)? strFolder + (lenStrFolder - 1) : strSep;
	strcat(strBuf, pSep);

	char *pFileName = strBuf + strlen(strBuf);
	strcpy(pFileName, arFldInFileName);

	//Reading field data from file:
	SRWLMagFld3D magFld3D;
	magFld3D.arBx = magFld3D.arBy = magFld3D.arBz = 0;
	magFld3D.arX = magFld3D.arY = magFld3D.arZ = 0; //if these pointers are not zero, regular mesh data will be ignored

	if(AuxReadInMagFld3D(&magFld3D, strBuf))
		cout << "Error reading 3D magnetic field data from file" << endl;
	cout << "done" << endl; 
	magFld3D.nRep = 1; //Entire ID

	SRWLMagFldC magFldCnt; //Container
	void *vArMagFld[] = {(void*)(&magFld3D)};
	magFldCnt.arMagFld = vArMagFld; //array of pointers to magnetic field elements
	magFldCnt.arMagFldTypes = "a"; //types of magnetic field elements in arMagFld array
	double auxArXcID[] = {xcID};
	double auxArYcID[] = {ycID};
	double auxArZcID[] = {zcID};
	magFldCnt.arXc = auxArXcID; //horizontal center positions of magnetic field elements in arMagFld array
	magFldCnt.arYc = auxArYcID; //vertical center positions of magnetic field elements in arMagFld array
	magFldCnt.arZc = auxArZcID; //longitudinal center positions of magnetic field elements in arMagFld array
	magFldCnt.nElem = 1; //number of magnetic field elements in arMagFld array

	//**********************Trajectory and Electron Beam structures
	SRWLPrtTrj partTraj;
	partTraj.partInitCond.x = 0.; //Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
	partTraj.partInitCond.y = 0.;
	partTraj.partInitCond.z = -129.027; //Initial Longitudinal Coordinate (set before the ID)
	partTraj.partInitCond.xp = 0.; //Initial Transverse Velocities
	partTraj.partInitCond.yp = 0.;
	partTraj.partInitCond.gamma = 17.5/0.51099890221e-03; //Relative Energy
	partTraj.partInitCond.relE0 = 1; //Electron Rest Mass
	partTraj.partInitCond.nq = -1; //Electron Charge

	partTraj.ctStart = 0.; //Start Time for the calculation
	partTraj.ctEnd = 270.;  //End Time
	partTraj.np = 537001; //Number of Points for Trajectory calculation

	//Allocating Trajectory arrays
	long totNumTrajPoints = partTraj.np*6;
	double *arTrajData = new double[totNumTrajPoints]; //one big array instead of 6 smaller ones
	double *t_arTrajData = arTrajData;
	//Setting pointers
	partTraj.arX = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arXp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arY = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arYp = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZ = t_arTrajData; t_arTrajData += partTraj.np;
	partTraj.arZp = t_arTrajData;

	//SRWLIB function call (to run calculations)
	cout << "   Performing trajectory calculation ... "; 
	ProcRes(srwlCalcPartTraj(&partTraj, &magFldCnt));
	cout << "done" << endl; 

	//Preparing full path to output ASCII file and Saving Trajectory data into it:
	strcpy(pFileName, strTrajOutFileName);
	cout << "   Saving trajectory data to a file ... "; 
	if(AuxSaveTrajData(&partTraj, strBuf)) 
		cout << "Error saving trajectory data to file" << endl;
	cout << "done" << endl; 

	//Electron Beam
	SRWLPartBeam elecBeam;
	elecBeam.Iavg = 0.5; //Average Current [A]
	elecBeam.partStatMom1 = partTraj.partInitCond;

	//**********************Precision parameters for SR calculation
	int meth = 1; //SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
	double relPrec = 0.01; //relative precision
	double zStartInteg = 0; //-129.029 #part.z - 0.1 #longitudinal position to start integration (effective if < zEndInteg)
	double zEndInteg = 0; //129.029 #part.z + 5.3 #longitudinal position to finish integration (effective if > zStartInteg)
	double sampFactNxNyForProp = 0; //sampling factor for adjusting nx, ny (effective if > 0)
	double arPrecPar[] = {meth, relPrec, zStartInteg, zEndInteg, partTraj.np, 0, sampFactNxNyForProp};

	//**********************Wavefront
	SRWLWfr wfr1; //For spectrum vs photon energy
	wfr1.mesh.ne = 5000; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr1.mesh.nx = wfr1.mesh.ny = 1;
	wfr1.mesh.zStart = 300.; //Longitudinal Position [m] at which SR has to be calculated
	wfr1.mesh.eStart = 5700.; //Initial Photon Energy [eV]
	wfr1.mesh.eFin = 6100.; //Final Photon Energy [eV]
	wfr1.mesh.xStart = 0.; //Initial Horizontal Position [m]
	wfr1.mesh.xFin = 0.; //Final Horizontal Position [m]
	wfr1.mesh.yStart = 0.; //Initial Vertical Position [m]
	wfr1.mesh.yFin = 0.; //Final Vertical Position [m]
	wfr1.partBeam = elecBeam;
	wfr1.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr1.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	long numTot = wfr1.mesh.ne*wfr1.mesh.nx*wfr1.mesh.ny*2;
	wfr1.arEx = (char*)(new float[numTot]); //horizontal and vertical electric field component arrays
	wfr1.arEy = (char*)(new float[numTot]);
	wfr1.arElecPropMatr = new double[20];
	wfr1.arWfrAuxData = new double[30];
	wfr1.arMomX = new double[11*wfr1.mesh.ne];
	wfr1.arMomY = new double[11*wfr1.mesh.ne];

	SRWLWfr wfr2; //For intensity distribution at fixed photon energy
	wfr2.mesh.ne = 1; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr2.mesh.nx = wfr2.mesh.ny = 61;
	wfr2.mesh.zStart = 300.; //Longitudinal Position [m] at which SR has to be calculated
	wfr2.mesh.eStart = 5915.5; //Initial Photon Energy [eV]
	wfr2.mesh.eFin = 5915.5; //Final Photon Energy [eV]
	wfr2.mesh.xStart = -0.0006; //Initial Horizontal Position [m]
	wfr2.mesh.xFin = 0.0006; //Final Horizontal Position [m]
	wfr2.mesh.yStart = -0.0006; //Initial Vertical Position [m]
	wfr2.mesh.yFin = 0.0006; //Final Vertical Position [m]
	wfr2.partBeam = elecBeam;
	wfr2.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr2.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	numTot = wfr2.mesh.ne*wfr2.mesh.nx*wfr2.mesh.ny*2;
	wfr2.arEx = (char*)(new float[numTot]); //horizontal and vertical electric field component arrays
	wfr2.arEy = (char*)(new float[numTot]);
	wfr2.arElecPropMatr = new double[20];
	wfr2.arWfrAuxData = new double[30];
	wfr2.arMomX = new double[11];
	wfr2.arMomY = new double[11];

	//**********************Calculation (SRWLIB function calls) and post-processing
	cout << "   Performing SR Spectrum vs Photon Energy calculation ... ";
	ProcRes(srwlCalcElecFieldSR(&wfr1, 0, &magFldCnt, arPrecPar));
	cout << "done" << endl;
	cout << "   Extracting Intensity from calculated Electric Field ... ";
	float *arI1 = new float[wfr1.mesh.ne]; //array to take spectrum data
	ProcRes(srwlCalcIntFromElecField((char*)arI1, &wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart));
	cout << "done" << endl;

	cout << "   Performing calculation of Electric Field at different positions ... ";
	ProcRes(srwlCalcElecFieldSR(&wfr2, 0, &magFldCnt, arPrecPar));
	cout << "done" << endl;
	cout << "   Extracting Intensity from calculated Electric Field ... ";
	float *arI2 = new float[wfr2.mesh.nx*wfr2.mesh.ny]; //"flat" array to take 2D intensity data
	ProcRes(srwlCalcIntFromElecField((char*)arI2, &wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0));
	cout << "done" << endl;

	//**********************Preparing full path to output ASCII file and Saving Intensity data into it:
	strcpy(pFileName, strIntOutFileName1);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI1, wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne, wfr1.mesh.xStart, wfr1.mesh.xFin, wfr1.mesh.nx, wfr1.mesh.yStart, wfr1.mesh.yFin, wfr1.mesh.ny, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl;

	strcpy(pFileName, strIntOutFileName2);
	cout << "   Saving intensity data to files ... ";
	if(AuxSaveIntensData(arI2, wfr2.mesh.eStart, wfr2.mesh.eFin, wfr2.mesh.ne, wfr2.mesh.xStart, wfr2.mesh.xFin, wfr2.mesh.nx, wfr2.mesh.yStart, wfr2.mesh.yFin, wfr2.mesh.ny, strBuf))
		cout << "Error saving resulting intensity data to file" << endl;
	cout << "done" << endl; 

	//**********************Deallocating memory
	delete[] magFld3D.arBx;
	delete[] magFld3D.arBy;
	delete[] magFld3D.arBz;
	delete[] arTrajData;
	delete[] wfr1.arEx;
	delete[] wfr1.arEy;
	delete[] wfr1.arWfrAuxData;
	delete[] wfr1.arMomX;
	delete[] wfr1.arMomY;
	delete[] wfr2.arEx;
	delete[] wfr2.arEy;
	delete[] wfr2.arWfrAuxData;
	delete[] wfr2.arMomX;
	delete[] wfr2.arMomY;
	return 0;
}

/************************************************************************//**
 * Main: illustrates and tests the basic functionality of SRW Library
 ***************************************************************************/
int main(int argc, char* argv[])
{
	char sVersNoSRW[1024], sVersNoSRWLIB[1024];
	ProcRes(srwlUtiVerNo(sVersNoSRW, 1));
	ProcRes(srwlUtiVerNo(sVersNoSRWLIB, 2));

	cout << "SRW Version: " << sVersNoSRW << endl; 
	cout << "SRWLIB Version: " << sVersNoSRWLIB << endl; 

	int numExamplesAvailable = 5;
	char *strPathProgr = argv[0];
	char *strPathExamples = 0;
	int exampleNo = -1; //run all examples by default
	char strBuf[2048];
	strBuf[0] = '\0';

	//Setting "callback" function pointer
	srwlUtiSetWfrModifFunc(&ModifySRWLWfr);

	if(argc > 1) exampleNo = atoi(argv[1]);
	if(argc > 2) strPathExamples = argv[2];
	else
	{
		char *pLastSep = strrchr(strPathProgr, '/');
		if(pLastSep == 0) pLastSep = strrchr(strPathProgr, '\\');

		strPathExamples = strBuf;
		if(pLastSep != 0) 
		{
			int nChar = (int)(pLastSep - strPathProgr + 1);
			strncpy(strBuf, strPathProgr, nChar);
			strBuf[nChar] = '\0';
		}
	}

	if(strPathExamples != 0)
	{
		if(exampleNo < 0)
		{//run all examples
			if(SRWLIB_Example01(strPathExamples)) cout << "Example # 1 was not executed correctly" << endl; 
			if(SRWLIB_Example02(strPathExamples)) cout << "Example # 2 was not executed correctly" << endl;
			if(SRWLIB_Example03(strPathExamples)) cout << "Example # 3 was not executed correctly" << endl;
			if(SRWLIB_Example04(strPathExamples)) cout << "Example # 4 was not executed correctly" << endl;
			if(SRWLIB_Example05(strPathExamples)) cout << "Example # 5 was not executed correctly" << endl;

			//to add new examples here
			//...
		}
		else if((exampleNo == 0) || (exampleNo > numExamplesAvailable))	cout << "Example # " << exampleNo << " doesn't exist" << endl; 
		else
		{
			//for(int iii=0; iii<200; iii++)
			//{
			bool failed = false;
			switch(exampleNo) {
				case 1:
					if(SRWLIB_Example01(strPathExamples)) failed = true; break;
				case 2:
					if(SRWLIB_Example02(strPathExamples)) failed = true; break;
				case 3:
					if(SRWLIB_Example03(strPathExamples)) failed = true; break;
				case 4:
					if(SRWLIB_Example04(strPathExamples)) failed = true; break;
				case 5:
					if(SRWLIB_Example05(strPathExamples)) failed = true; break;
				//to add new examples here
			}
			if(failed) cout << "Example # " << exampleNo << " was not executed correctly" << endl;
			//}
		}
	}
	else cout << "Path to SRWLIB Examples was not specified" << endl;

	cout << "Press Enter to exit" << endl;
	getchar();
	return 0;
}

/***************************************************************************/
