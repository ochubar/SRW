
//+++++++++++++++++++++++++++++++++++++++
//
//Test calculation of Fresnel Diffraction on a Slit
//Usage: 
//wTestFresnDif = cabs(srwWfrPointSourceFresnelDif(90000, x, 0, 0, 0, 10e-06, 50e-06, 4.365, 5.703))^2
//+++++++++++++++++++++++++++++++++++++++
function/C srwWfrPointSourceFresnelDif(phEn, x, y, x0, y0, dx, dy, R, L)
variable phEn, x, y, x0, y0, dx, dy, R, L

variable halfWaveNumber = phEn*2.53384080189E+06
variable pp = halfWaveNumber*(1/R + 1/L)
variable sqrFact = sqrt(pp/(2*Pi))
variable invRpL = 1/(R + L)
variable xc = (x0*L + x*R)*invRpL, yc = (y0*L + y*R)*invRpL
variable argX1 = sqrFact*(dx - 2*xc), argX2 = sqrFact*(dx + 2*xc)
variable argY1 = sqrFact*(dy - 2*yc), argY2 = sqrFact*(dy + 2*yc)

variable Cx1 = fresnelCos(argX1), Cx2 = fresnelCos(argX2)
variable Sx1 = fresnelSin(argX1), Sx2 = fresnelSin(argX2)
variable Cy1 = fresnelCos(argY1), Cy2 = fresnelCos(argY2)
variable Sy1 = fresnelSin(argY1), Sy2 = fresnelSin(argY2)

variable/C Fx = cmplx(Cx1 + Cx2, Sx1 + Sx2)
variable/C Fy = cmplx(Cy1 + Cy2, Sy1 + Sy2)
return Fx*Fy
end

//+++++++++++++++++++++++++++++++++++++++
//
//(???) Test calculation of Fresnel Diffraction on a Slit
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropagTestFresnelDif(StoName, ElecName, Wfr, BL, PrecParM, MaxPrt)
string StoName=srwUtiGetValS("StoName", "Stk", "SrwWfrPropagStokesMultiE")
string ElecName=SrwElecName+SrwElecType
string Wfr=SrwRadName+SrwRadType
string BL=SrwBliLast+SrwBeamlineType
variable PrecParM=srwUtiGetValN("PrecParM", 1, "SrwWfrPropagStokesMultiE")
variable MaxPrt=srwUtiGetValN("MaxPrt", 10000, "SrwWfrPropagStokesMultiE") 
prompt StoName, "Name of the Stokes structure"
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt Wfr, "Wavefront",popup Wavelist("*"+SrwRadType ,";", "")
prompt BL, "Optical component",popup Wavelist("*"+SrwBeamlineType ,";", "")
prompt PrecParM, "Multi-e propag. precision"
prompt MaxPrt, "Max. number of macro-particles"
Silent 1						|	Propagating the Wavefront ...
PauseUpdate

SrwBliLast = BL[0,strlen(BL)-strlen(SrwBeamlineType)-1]
SrwElecName = ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwRadName = Wfr[0,strlen(Wfr)-strlen(SrwRadType)-1]
srwUtiSetValS("StoName", StoName, "SrwWfrPropagStokesMultiE")
srwUtiSetValN("PrecParM", PrecParM, "SrwWfrPropagStokesMultiE")
srwUtiSetValN("MaxPrt", MaxPrt, "SrwWfrPropagStokesMultiE")

String AuxObsName = "AuxObs"
String AuxObs = AuxObsName + SrwSmpType
SrwSmpCreate(AuxObsName, 1.)
SrwSmpScanXZE(AuxObs, 0, 1, 10, 0, 1, 10, 19., 19., 1)
SrwStoPrepSimple(AuxObs,StoName,2)

Make/D/O/N=7 waveprec
waveprec[0]=PrecParM  // Rel. Single-E Propag. Prec.
waveprec[1]=MaxPrt

srRadPropagStokesMultiE($ElecName, $Wfr, $BL, waveprec, $(StoName + SrwStoType))

KillWaves/Z  waveprec, AuxObs
end
