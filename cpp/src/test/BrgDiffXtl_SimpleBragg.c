// John P. Sutter
// 26 July 2011
// This program calculates the final wave vector and final electric field
// amplitudes for a given incident wave vector with initial electric field
// amplitudes. Dynamical diffraction theory for a perfect crystal is applied.
// Only the 2-beam Bragg case with diffraction plane normal to the crystal
// surface is considered, but the user may view the transmitted beam as well
// as the reflected beam. Equations are from W. H. Zachariasen, "Theory of
// X-ray Diffraction in Crystals" Chapter III.

#include<iostream>
#include<complex>
#include<limits>
using namespace std;

int main () {

// LIST OF USEFUL CONSTANTS
 
// conversion factor: structure factor -> electric susceptibility Fourier component
     const double epmcA = 8.969785162E-06; 

// conversion factor: wavelength <-> energy
     const double eAconv = 12398.4193009;

// Silicon lattice parameters (in Angstroms):
//  300 K
     const double aSiART = 5.43108124; 
//  120 K: minimum value
     const double aSiALT = 5.42965775;

// Germanium lattice parameters for 70Ge (in Angstroms):
//  300 K
     const double aGeART = 5.6579830;
//  75 K
     const double aGeALT = 5.6525930;

// Diamond lattice parameters (in Angstroms):
//  298 K
     const double aDiART = 3.56712;
//  77 K
     const double aDiALT = 3.56681;

// Sapphire lattice parameters for hexagonal axes (in Angstroms):
// 300 K
     const double aSaART = 4.7593240;
     const double cSaART = 12.9919261;
// 80 K
     const double aSaALT = 4.7563257;
     const double cSaALT = 12.9820703;

// Quartz lattice parameters for hexagonal axes (in Angstroms):
// 296 K
     const double aQuART = 4.9141;
     const double cQuART = 5.4060;
// 78 K
     const double aQuALT = 4.9030;
     const double cQuALT = 5.3999; 

// END OF CONSTANTS

// USER INPUTS

// Input #1a: crystal material
     int xtltyp = 1;
// Input #1b: choice of room-temperature or low-temperature lattice parameters:
// rmtmp = 0 for low temperature, 1 for room temperature
     int rmtmp = 1;

// Input #2: photon energy in eV
     double EceV = 11183.9195;

// Input #3: Miller indices of Bragg reflection. See Kittel for definitions.
     int hMilND = 4;
     int kMilND = 4;
     int lMilND = 4;

// Input #4: Asymmetry angle in diffraction plane: (+) for grazing incidence and
// (-) for grazing exit. In degrees.
     double alphdg = 0.0;

// Input #5: Input mode type. Two modes are available:
// AUTO: User inputs the deviation from the Bragg angle, which is calculated
//       automatically. (ptMode = 0)
// ABSOLUTE: User inputs the absolute angle relative to the diffracting planes.
//       This is not equal to the pitch rotation angle of the crystal if the
//       asymmetry angle is nonzero. (ptMode = 1)
     int ptMode = 0; 
// Input #5a: energy at which crystal's Bragg condition should be fulfilled when
//   the deviation (input #5b) is zero. This need not equal the central energy
//   of the wavefront! Given in eV.
     double EadjeV = 11183.9195;
// Input #5b: deviation from kinematic Bragg angle at E given by input #5a, in urad
     double dthur = 0.0;
// Input #5c: absolute incidence angle from diffracting planes in degrees. This is
// used only if ptMode = 1.
     double thdg = 42.;

// Input #6: roll angle in degrees
     double chidg = 0.0;

// Input #7: yaw angle in degrees
     double phidg = 0.0;

// Input #8: crystal thickness in microns
     double thicum = 500.;

// Input #9: Whether to calculate the transmitted beam as well as the diffracted
// beam. itrans = 0 for NO, itrans = 1 for YES.
     int itrans = 1;

// Input #10: structure factor for k=0. NOTE: This is hard-coded
// in this version, but should ideally be calculated from the atomic positions
// in the lattice.
     complex<double> f0c(113.3969,1.3726);

// Input #11: structure factor for k=H (diffracted beam). The note for Input
// #10 also applies here.
     complex<double> fhc(-34.4013,-1.1446);

// Input #12: structure factor for k=-H. The note for Input #10 also applies
// here.
     complex<double> fmhc(-34.4013,-1.1446);

// END OF USER INPUTS

// PRELIMINARY CALCULATIONS FROM USER INPUTS
     double aSiA,aGeA,aDiA,aSaA,cSaA,aQuA,cQuA,volA3;
     double alamA,k0Ai,dA,alphrd,HXAi[3],thptrd,chird,phird;
     double aladjA,thBrd,thrd;
     double pi,DBLEMax,logDBMax;
     complex<double> iC(0.,1.),psi0c,psihc,psimhc;

     pi = 4.*atan(1.);
     DBLEMax = numeric_limits<double>::max();
     logDBMax = log(DBLEMax);

// From Input #1a,b: volume of crystal lattice unit cell
     if (xtltyp == 1){ // silicon
        if (rmtmp == 0) // low temperature
           aSiA = aSiALT;
        else // high temperature
           aSiA = aSiART;
        volA3 = pow(aSiA,3); 
     }
     else if (xtltyp == 2){ // germanium
        if (rmtmp == 0)
           aGeA = aGeALT;
        else
           aGeA = aGeART;
        volA3 = pow(aGeA,3);
     }
     else if (xtltyp == 3){ // diamond
        if (rmtmp == 0)
           aDiA = aDiALT;
        else
           aDiA = aDiART;
        volA3 = pow(aDiA,3);
     }
     else if (xtltyp == 4){ // sapphire
        if (rmtmp == 0){
           aSaA = aSaALT;
           cSaA = cSaALT;
        }
        else{
           aSaA = aSaART;
           cSaA = cSaART; 
        }
        volA3 = 0.5*sqrt(3.0)*aSaA*aSaA*cSaA;
     }
     else if (xtltyp == 5){ // quartz
        if (rmtmp = 0){
           aQuA = aQuALT;
           cQuA = cQuALT;
        }
        else{
           aQuA = aQuART;
           cQuA = cQuART;
        }
        volA3 = 0.5*sqrt(3.0)*aQuA*aQuA*cQuA;
     }

// From Input #2: conversion of photon energy to wavelength, and conversion
// of structure factor for k=0 to Fourier component of electric susceptibility
     alamA = eAconv/EceV;
     k0Ai = 1./alamA;
     psi0c = -epmcA*alamA*alamA*f0c/volA3;    

// From Input #3: interplanar spacing, and conversion of structure factor for
// k=H and k=-H to corresponding Fourier components of electric susceptibility
     if (xtltyp == 1) // silicon
        dA = aSiA/sqrt( hMilND*hMilND + kMilND*kMilND + lMilND*lMilND );
     else if (xtltyp == 2) // germanium
        dA = aGeA/sqrt( hMilND*hMilND + kMilND*kMilND + lMilND*lMilND );
     else if (xtltyp == 3) // diamond
        dA = aDiA/sqrt( hMilND*hMilND + kMilND*kMilND + lMilND*lMilND );
     else if (xtltyp == 4) // sapphire
        dA = 1./( sqrt( 4.*( hMilND*hMilND + kMilND*kMilND + hMilND*kMilND )/
             (3.*aSaA*aSaA) + (lMilND*lMilND)/(cSaA*cSaA) ) );
     else if (xtltyp == 5) // quartz 
        dA = 1./( sqrt( 4.*( hMilND*hMilND + kMilND*kMilND + hMilND*kMilND )/
             (3.*aQuA*aQuA) + (lMilND*lMilND)/(cQuA*cQuA) ) );
     psihc  = -epmcA*alamA*alamA*fhc/volA3;    
     psimhc = -epmcA*alamA*alamA*fmhc/volA3;    
     
// From Input #4: asymmetry angle in radians (for trigonometric functions)
// and length of reciprocal lattice vector of Bragg reflection
     alphrd = alphdg*pi/180.;
     HXAi[0] =  0.;
     HXAi[1] =  cos(alphrd)/dA;
     HXAi[2] = -sin(alphrd)/dA;

// From Input #5: the photon energy used to orient the crystal, the Bragg
// angle of the reflection, and the pitch rotation angle of the crystal 
// in radians. This is the angle of rotation of the crystal coordinate 
// system relative to the lab frame.
     if (ptMode == 0){
        aladjA = eAconv/EadjeV;
        thBrd = asin(aladjA/(2.*dA)); 
        thrd = thBrd + dthur*1.E-06;
     }
     else
        thrd = thdg*pi/180.;
     thptrd = thrd - alphrd;
     
// From Input #6: conversion of roll angle to radians
     chird = chidg*pi/180.;

// From Input #7: conversion of yaw angle to radians
     phird = phidg*pi/180.;

// END PRELIMINARY CALCULATIONS FROM USER INPUTS

// PRINT STATEMENTS FOR TESTS: MAY BE REMOVED
     printf("pi = %f\n",pi);
     printf("iC = %f,%f\n",real(iC),imag(iC));
     printf("DBLEMax = %e\n",DBLEMax);
     printf("logDBMax = %e\n",logDBMax);

     printf("Real and imaginary parts of psi0c:\n");
     printf("%e,%e\n",real(psi0c),imag(psi0c));
     
     printf("Real and imaginary parts of psihc:\n");
     printf("%e,%e\n",real(psihc),imag(psihc));
  
     printf("Real and imaginary parts of psimhc:\n");
     printf("%e,%e\n",real(psimhc),imag(psimhc));

     printf("Interplanar spacing = %f\n",dA);
     printf("Recip. vector = %f,%f,%f\n",HXAi[0],HXAi[1],HXAi[2]);

     printf("Bragg angle = %f\n",thBrd*180./pi);
     printf("Incidence angle = %f\n",thptrd*180./pi);
// END PRINT STATEMENTS FOR TESTS

// TRANSFORMATION MATRIX
// RXtLab: 3x3 orthogonal matrix that converts components of a 3x1 vector
// in crystal coordinates to components in lab coordinates.
// RLabXt: transpose of RXtLab: converts components of a 3x1 vector in lab
// coordinates to components in crystal coordinates.

     int i,j;
     double RXtLab[3][3],RLabXt[3][3];

     RXtLab[0][0] =  cos(chird)*cos(phird);
     RXtLab[0][1] = -sin(chird); 
     RXtLab[0][2] =  cos(chird)*sin(phird);
     RXtLab[1][0] =  cos(thptrd)*sin(chird)*cos(phird)-sin(thptrd)*sin(phird);
     RXtLab[1][1] =  cos(thptrd)*cos(chird);
     RXtLab[1][2] =  cos(thptrd)*sin(chird)*sin(phird)+sin(thptrd)*cos(phird);
     RXtLab[2][0] = -sin(thptrd)*sin(chird)*cos(phird)-cos(thptrd)*sin(phird);
     RXtLab[2][1] = -sin(thptrd)*cos(chird);
     RXtLab[2][2] = -sin(thptrd)*sin(chird)*sin(phird)+cos(thptrd)*cos(phird);

     for (i = 0; i < 3; i++){
         for (j = 0; j < 3; j++)
             RLabXt[i][j] = RXtLab[j][i]; 
     }

// END OF TRANSFORMATION MATRIX

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("RXtLab =\n");
     printf("%e  %e  %e\n",RXtLab[0][0],RXtLab[0][1],RXtLab[0][2]);
     printf("%e  %e  %e\n",RXtLab[1][0],RXtLab[1][1],RXtLab[1][2]);
     printf("%e  %e  %e\n",RXtLab[2][0],RXtLab[2][1],RXtLab[2][2]);
     printf("RLabXt =\n");
     printf("%e  %e  %e\n",RLabXt[0][0],RLabXt[0][1],RLabXt[0][2]);
     printf("%e  %e  %e\n",RLabXt[1][0],RLabXt[1][1],RLabXt[1][2]);
     printf("%e  %e  %e\n",RLabXt[2][0],RLabXt[2][1],RLabXt[2][2]);
// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED

// DYNAMICAL DIFFRACTION CALCULATIONS
// NOTE: This is a loop over various values of the wavevector components
// (kx,ky).
     
     int ik,nk = 1001;
     double e1[3],e2[3],e1X[3],e2X[3],PolTrn[2][2];
     double kc0XAi[3],ucn0X[3],sg0X[3],pi0X[3],sgmag;
     double kcHtXA[3],kcHXAi[3],kcHt2;
     double ucnHX[3],sgHX[3],piHX[3],sgH[3],piH[3];
     double kxAi,kyAi,kzAi,k0wAi[3],k0wXAi[3],u0[3],u0X[3];
     double kHtXAi[3],kHt2,kmHXAi[3],kmHAi[3];
     double gamma0,bee,dotkH,Adev,cos2t;
     double dangur,uref[3],urkmH,kangur,ukcrAi[3],mgkmk0;
     double mgDHsg,mgDHpi,mgEHs,mgEHp,mgE0ts,mgE0tp;
     complex<double> EIn12s,EIn12p,EInSPs,EInSPp,zeeC,queC,sqrqzC;
     complex<double> del1C,del2C,x1C,x2C,ph1C,ph2C,Cph1C,Cph2C;
     complex<double> DHsgC,DHpiC,D0trsC,D0trpC;
     complex<double> EHSPCs,EHSPCp,E0tSPs,E0tSPp;
     FILE *TestK, *TestRef, *TestAmp, *TestTrns;
// THESE FILES ARE FOR TESTS ONLY AND MAY BE CHANGED IN THE FINAL VERSION!
     TestK = fopen("TestK.out","w");
     TestRef = fopen("TestRef.out","w");
     TestAmp = fopen("TestAmp.out","w");
     TestTrns = fopen("TestTrns.out","w");
// Polarization vectors e1,e2 in lab frame:
// z along beam, y upward, x = y x z
// xz-plane is horizontal, yz-plane is vertical
// NOTE: We neglect dependence of polarization vectors on (kx,ky).
// If the angular range of the wave vectors is under 1 mrad, this should
// be OK.
     e1[0] = 1.;
     e1[1] = 0.;
     e1[2] = 0.; 
     e2[0] = 0.;
     e2[1] = 1.;
     e2[2] = 0.;
// Conversion of polarization vector components to crystal frame:
     e1X[0] = RLabXt[0][0];
     e1X[1] = RLabXt[1][0];
     e1X[2] = RLabXt[2][0];
     e2X[0] = RLabXt[0][1];
     e2X[1] = RLabXt[1][1];
     e2X[2] = RLabXt[2][1]; 

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("e1X = \n");
     printf("%e\n",e1X[0]);
     printf("%e\n",e1X[1]);
     printf("%e\n",e1X[2]);
     printf("e2X = \n");
     printf("%e\n",e2X[0]);
     printf("%e\n",e2X[1]);
     printf("%e\n",e2X[2]);
// END PRINT STATEMENTS FOR TESTS

// Calculation of new polarization vectors in the crystal frame, with the
// sigma vector parallel to the diffracting atomic planes of the crystal.
// In this simple case, the sigma vector will also be parallel to the 
// crystal surface.
// This is an especially convenient choice because the sigma and pi components
// of the wavefields inside the crystal can be treated separately.
// Again, we neglect the dependence on (kx,ky).

     kc0XAi[0] = k0Ai*RLabXt[0][2];
     kc0XAi[1] = k0Ai*RLabXt[1][2];
     kc0XAi[2] = k0Ai*RLabXt[2][2]; 
     ucn0X[0] = RLabXt[0][2];
     ucn0X[1] = RLabXt[1][2];
     ucn0X[2] = RLabXt[2][2];
     sg0X[0] = HXAi[1]*ucn0X[2] - HXAi[2]*ucn0X[1];
     sg0X[1] = HXAi[2]*ucn0X[0] - HXAi[0]*ucn0X[2];
     sg0X[2] = HXAi[0]*ucn0X[1] - HXAi[1]*ucn0X[0];
     sgmag = sqrt(sg0X[0]*sg0X[0]+sg0X[1]*sg0X[1]+sg0X[2]*sg0X[2]);
     sg0X[0] = sg0X[0]/sgmag;
     sg0X[1] = sg0X[1]/sgmag;
     sg0X[2] = sg0X[2]/sgmag;
     pi0X[0] = ucn0X[1]*sg0X[2] - ucn0X[2]*sg0X[1];
     pi0X[1] = ucn0X[2]*sg0X[0] - ucn0X[0]*sg0X[2];
     pi0X[2] = ucn0X[0]*sg0X[1] - ucn0X[1]*sg0X[0];

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("ucn0X =\n");
     printf("%e\n",ucn0X[0]);
     printf("%e\n",ucn0X[1]);
     printf("%e\n",ucn0X[2]);
     printf("sg0X =\n");
     printf("%e\n",sg0X[0]);
     printf("%e\n",sg0X[1]);
     printf("%e\n",sg0X[2]);
     printf("pi0X =\n");
     printf("%e\n",pi0X[0]);
     printf("%e\n",pi0X[1]);
     printf("%e\n",pi0X[2]);
// END PRINT STATEMENTS FOR TESTS

// Calculate the 2x2 transformation matrix of the polarizations from 
// (e1X,e2X) to (sg0X,pi0X).

     PolTrn[0][0] = e1X[0]*sg0X[0] + e1X[1]*sg0X[1] + e1X[2]*sg0X[2];
     PolTrn[0][1] = e2X[0]*sg0X[0] + e2X[1]*sg0X[1] + e2X[2]*sg0X[2];
     PolTrn[1][0] = e1X[0]*pi0X[0] + e1X[1]*pi0X[1] + e1X[2]*pi0X[2];
     PolTrn[1][1] = e2X[0]*pi0X[0] + e2X[1]*pi0X[1] + e2X[2]*pi0X[2];

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("PolTrn =\n");
     printf("%e  %e\n",PolTrn[0][0],PolTrn[0][1]);
     printf("%e  %e\n",PolTrn[1][0],PolTrn[1][1]);
// END OF PRINT STATEMENTS FOR TESTS

// Calculate polarization vectors in the crystal frame for the diffracted
// beam, ignoring the dependence on (kx,ky):

     kcHtXA[0] = kc0XAi[0] + HXAi[0];
     kcHtXA[1] = 0.;
     kcHtXA[2] = kc0XAi[2] + HXAi[2]; 
     kcHt2 = kcHtXA[0]*kcHtXA[0] + kcHtXA[2]*kcHtXA[2];
     kcHXAi[0] = kcHtXA[0];
     kcHXAi[1] = sqrt(k0Ai*k0Ai-kcHt2);
     kcHXAi[2] = kcHtXA[2];
     ucnHX[0] = kcHXAi[0]/k0Ai;
     ucnHX[1] = kcHXAi[1]/k0Ai;
     ucnHX[2] = kcHXAi[2]/k0Ai;
     sgHX[0] = sg0X[0];
     sgHX[1] = sg0X[1];
     sgHX[2] = sg0X[2];
     piHX[0] = ucnHX[1]*sgHX[2] - ucnHX[2]*sgHX[1];
     piHX[1] = ucnHX[2]*sgHX[0] - ucnHX[0]*sgHX[2];
     piHX[2] = ucnHX[0]*sgHX[1] - ucnHX[1]*sgHX[0];

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("ucnHX =\n");
     printf("%e\n",ucnHX[0]);
     printf("%e\n",ucnHX[1]);
     printf("%e\n",ucnHX[2]);
     printf("sgHX =\n");
     printf("%e\n",sgHX[0]);
     printf("%e\n",sgHX[1]);
     printf("%e\n",sgHX[2]);
     printf("piHX =\n");
     printf("%e\n",piHX[0]);
     printf("%e\n",piHX[1]);
     printf("%e\n",piHX[2]);
// END PRINT STATEMENTS FOR TESTS
 
// Transform the diffracted beam polarization vectors back into the lab frame:

     sgH[0] = RXtLab[0][0]*sgHX[0] + RXtLab[0][1]*sgHX[1] + RXtLab[0][2]*sgHX[2];
     sgH[1] = RXtLab[1][0]*sgHX[0] + RXtLab[1][1]*sgHX[1] + RXtLab[1][2]*sgHX[2];
     sgH[2] = RXtLab[2][0]*sgHX[0] + RXtLab[2][1]*sgHX[1] + RXtLab[2][2]*sgHX[2];
     piH[0] = RXtLab[0][0]*piHX[0] + RXtLab[0][1]*piHX[1] + RXtLab[0][2]*piHX[2];
     piH[1] = RXtLab[1][0]*piHX[0] + RXtLab[1][1]*piHX[1] + RXtLab[1][2]*piHX[2];
     piH[2] = RXtLab[2][0]*piHX[0] + RXtLab[2][1]*piHX[1] + RXtLab[2][2]*piHX[2];

// PRINT STATEMENTS FOR TESTS: CAN BE REMOVED
     printf("sgH =\n");
     printf("%e\n",sgH[0]);
     printf("%e\n",sgH[1]);
     printf("%e\n",sgH[2]);
     printf("piH =\n");
     printf("%e\n",piH[0]);
     printf("%e\n",piH[1]);
     printf("%e\n",piH[2]);
// END PRINT STATEMENTS FOR TESTS

// Begin loop over wave vectors:

     for (ik = 0; ik < nk; ik++){
// Wave vector components: hard-coded for this test only!
         kxAi = 0.0;
         kyAi = -0.00025 + 0.00050*(double)(ik)/(double)(nk-1);
// Electric field amplitudes in k-space: hard-coded for this test only!
         EIn12s = complex<double>(0.,0.);
         EIn12p = complex<double>(1.,0.);  
// k0wvAi = incident beam wave vector (kx,ky,kz).
         kzAi = sqrt(k0Ai*k0Ai - kxAi*kxAi - kyAi*kyAi);
         k0wAi[0] = kxAi;
         k0wAi[1] = kyAi;
         k0wAi[2] = kzAi;
         u0[0] = kxAi/k0Ai;
         u0[1] = kyAi/k0Ai;
         u0[2] = kzAi/k0Ai;
// Conversion of wave vector components to crystal frame:
         k0wXAi[0] = RLabXt[0][0]*k0wAi[0] + RLabXt[0][1]*k0wAi[1] + RLabXt[0][2]*k0wAi[2];
         k0wXAi[1] = RLabXt[1][0]*k0wAi[0] + RLabXt[1][1]*k0wAi[1] + RLabXt[1][2]*k0wAi[2];
         k0wXAi[2] = RLabXt[2][0]*k0wAi[0] + RLabXt[2][1]*k0wAi[1] + RLabXt[2][2]*k0wAi[2];
         u0X[0] = k0wXAi[0]/k0Ai;
         u0X[1] = k0wXAi[1]/k0Ai;
         u0X[2] = k0wXAi[2]/k0Ai;
// Convert components of incident beam polarization from (e1X,e2X) to (sg0X,pi0X):
         EInSPs = PolTrn[0][0]*EIn12s + PolTrn[0][1]*EIn12p;
         EInSPp = PolTrn[1][0]*EIn12s + PolTrn[1][1]*EIn12p; 
// Calculate the wave vector for the diffracted beam. Note that here the index of
// refraction corrections are included.
         kHtXAi[0] = k0wXAi[0] + HXAi[0];
         kHtXAi[1] = 0.;
         kHtXAi[2] = k0wXAi[2] + HXAi[2]; 
         kHt2 = kHtXAi[0]*kHtXAi[0] + kHtXAi[2]*kHtXAi[2];
         kmHXAi[0] = kHtXAi[0];
         kmHXAi[1] = sqrt(k0Ai*k0Ai-kHt2);
         kmHXAi[2] = kHtXAi[2];
// Calculate direction cosine gamma0,reflection asymmetry parameter bee,
// deviation parameter Adev and normalized deviation parameter zeeC (as
// defined by Zachariasen gamma0, b, alpha and z.)
         gamma0 = -u0X[1];
         bee = 1./(1. + HXAi[1]/k0wXAi[1]);
         dotkH = k0wXAi[0]*HXAi[0] + k0wXAi[1]*HXAi[1] + k0wXAi[2]*HXAi[2];
         Adev = ( 2.*dotkH + 1./(dA*dA) )/(k0Ai*k0Ai);
         zeeC = 0.5*( (1.-bee)*psi0c + bee*Adev );
// Calculate the complex reflectivity DHsgC for sigma polarization:
         queC = bee*psihc*psimhc;
         sqrqzC = sqrt(queC + zeeC*zeeC); 
         del1C = 0.5*( psi0c - zeeC + sqrqzC );
         del2C = 0.5*( psi0c - zeeC - sqrqzC );
         x1C = (-zeeC+sqrqzC)/psimhc;
         x2C = (-zeeC-sqrqzC)/psimhc;
         ph1C = -2.*pi*iC*k0Ai*del1C*(thicum*1.E+04)/gamma0;
         ph2C = -2.*pi*iC*k0Ai*del2C*(thicum*1.E+04)/gamma0;
         if (real(ph1C) > logDBMax)
            DHsgC = x2C;
         else if (real(ph2C) > logDBMax)
            DHsgC = x1C;
         else {
            Cph1C = exp(ph1C);
            Cph2C = exp(ph2C);
            DHsgC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C); 
         } 
         if (itrans == 1) {
//    calculate the complex reflectivity D0trsC of the transmitted beam.
            if (real(ph1C) > logDBMax){
               Cph2C = exp(ph2C);
               D0trsC = -Cph2C*(x2C-x1C)/x1C;
            }
            else if (real(ph2C) > logDBMax){
               Cph1C = exp(ph1C);
               D0trsC = +Cph1C*(x2C-x1C)/x2C;
            }
            else
               D0trsC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C);
         }
// Calculate the complex reflectivity DHpiC for pi polarization:
         cos2t = cos(2.*thrd);
         queC = bee*psihc*psimhc*cos2t*cos2t;
         sqrqzC = sqrt(queC + zeeC*zeeC); 
         del1C = 0.5*( psi0c - zeeC + sqrqzC );
         del2C = 0.5*( psi0c - zeeC - sqrqzC );
         x1C = (-zeeC+sqrqzC)/(psimhc*cos2t);
         x2C = (-zeeC-sqrqzC)/(psimhc*cos2t);
         ph1C = -2.*pi*iC*k0Ai*del1C*(thicum*1.E+04)/gamma0;
         ph2C = -2.*pi*iC*k0Ai*del2C*(thicum*1.E+04)/gamma0;
         if (real(ph1C) > logDBMax)
            DHpiC = x2C;
         else if (real(ph2C) > logDBMax)
            DHpiC = x1C;
         else {
            Cph1C = exp(ph1C);
            Cph2C = exp(ph2C);
            DHpiC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C); 
         } 
         if (itrans == 1) {
//    calculate the complex reflectivity D0trpC of the transmitted beam.
            if (real(ph1C) > logDBMax){
               Cph2C = exp(ph2C);
               D0trpC = -Cph2C*(x2C-x1C)/x1C;
            }
            else if (real(ph2C) > logDBMax){
               Cph1C = exp(ph1C);
               D0trpC = +Cph1C*(x2C-x1C)/x2C;
            }
            else
               D0trpC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C);
         }
// Calculate the diffracted amplitudes:
         EHSPCs = DHsgC*EInSPs;
         EHSPCp = DHpiC*EInSPp; 
         if (itrans == 1){
            E0tSPs = D0trsC*EInSPs;
            E0tSPp = D0trpC*EInSPp;
         }
// Convert the diffracted beam's wave vector from the crystal frame back to the lab frame.
         kmHAi[0] = RXtLab[0][0]*kmHXAi[0] + RXtLab[0][1]*kmHXAi[1] + RXtLab[0][2]*kmHXAi[2];
         kmHAi[1] = RXtLab[1][0]*kmHXAi[0] + RXtLab[1][1]*kmHXAi[1] + RXtLab[1][2]*kmHXAi[2];
         kmHAi[2] = RXtLab[2][0]*kmHXAi[0] + RXtLab[2][1]*kmHXAi[1] + RXtLab[2][2]*kmHXAi[2];
// For test output only:
         dangur = -asin(kyAi/k0Ai)*1.E+06;
         uref[0] = 0.;
         uref[1] = sin(2.*thBrd);
         uref[2] = cos(2.*thBrd);
         urkmH = uref[0]*kmHAi[0] + uref[1]*kmHAi[1] + uref[2]*kmHAi[2];
         if (urkmH > k0Ai)
            kangur = 0.;
         else
            kangur = acos(urkmH/k0Ai)*1.E+06;
         ukcrAi[0] = uref[1]*kmHAi[2] - uref[2]*kmHAi[1]; 
         ukcrAi[1] = uref[2]*kmHAi[0] - uref[0]*kmHAi[2]; 
         ukcrAi[2] = uref[0]*kmHAi[1] - uref[1]*kmHAi[0]; 
         if (ukcrAi[0] > 0.)
            kangur = -kangur;
         mgkmk0 = sqrt( kmHAi[0]*kmHAi[0] + kmHAi[1]*kmHAi[1] + kmHAi[2]*kmHAi[2] )/k0Ai;
         mgDHsg = abs(DHsgC);
         mgDHpi = abs(DHpiC);
         mgEHs = abs(EHSPCs);
         mgEHp = abs(EHSPCp);
         fprintf(TestK,"%e  %e  %e\n",dangur,mgkmk0,kangur);
         fprintf(TestRef,"%e  %e  %e  %e  %e\n",dangur,mgDHsg*mgDHsg,arg(DHsgC),mgDHpi*mgDHpi,arg(DHpiC));
         fprintf(TestAmp,"%e  %e  %e  %e  %e\n",dangur,mgEHs*mgEHs,arg(EHSPCs),mgEHp*mgEHp,arg(EHSPCp));
         if (itrans == 1){
            mgE0ts = abs(E0tSPs);
            mgE0tp = abs(E0tSPp);
            fprintf(TestTrns,"%e  %e  %e  %e  %e\n",dangur,mgE0ts*mgE0ts,arg(E0tSPs),mgE0tp*mgE0tp,arg(E0tSPp));
         }
     }

     fclose(TestK);
     fclose(TestRef);
     fclose(TestAmp);
     fclose(TestTrns);

// the return statement is required by some compilers though not all
     return 0;
}
