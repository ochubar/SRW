
//-----------------------------------------------
//Prepares special Phase Shift wave for focusing of arbitrary electric field wavefront
//-----------------------------------------------
proc SrwMiscPrepSpecZP(InE, OutTName, L, x0, z0)
String InE, OutTName
Variable L, x0, z0

Variable PhEn_eV = DimOffset($InE, 0)
Variable PhEn_keV = 0.001*PhEn_eV

Variable LambdaMic = 1.239854/PhEn_eV
Variable WaveNumbInvM = 2*Pi*(1.e+06)/LambdaMic

Variable nx = DimSize($InE, 1), nz = DimSize($InE, 2)

Variable xRange = DimDelta($InE, 1)*(nx - 1), zRange = DimDelta($InE, 2)*(nz - 1)
Variable xc = DimOffset($InE, 1) + 0.5*xRange, zc = DimOffset($InE, 2) + 0.5*zRange
Variable xc_mm = 1000.*xc, zc_mm = 1000.*zc, xRange_mm = 1000.*xRange, zRange_mm = 1000.*zRange

SrwOptThinGenTempl(OutTName,PhEn_keV,xc_mm,zc_mm,xRange_mm,zRange_mm,nx,nz,1)
OutTName += SrwOptThinGenWaveType

$OutTName = cmplx(1, imag(r2polar(srwMiscETrans($InE[0][p][q], WaveNumbInvM, L, x, y, x0, z0))))
end

//-----------------------------------------------
function/C srwMiscETrans(InE, WaveNumb, L, x, z, x0, z0)
Variable/C InE
Variable WaveNumb, L, x, z, x0, z0

Variable xr = x - x0, zr = z - z0
//Variable Ph = WaveNumb*sqrt(L*L + xr*xr + zr*zr)
Variable Ph = -WaveNumb*(xr*xr + zr*zr)/(2*L)

Variable/C ExpMult = cmplx(cos(Ph), sin(Ph))
if(InE != 0)
	return ExpMult/InE
else
	return ExpMult
endif
end

//-----------------------------------------------
//Prepares special Phase Shift wave for compensation of the wavefront differences from non-point coherent source
//-----------------------------------------------
proc SrwMiscPrepCompensZP(InE, OutTName, Lsrc, x0src, z0src, Ph0)
String InE, OutTName
Variable Lsrc, x0src, z0src, Ph0

Variable PhEn_eV = DimOffset($InE, 0)
Variable PhEn_keV = 0.001*PhEn_eV

Variable LambdaMic = 1.239854/PhEn_eV
Variable WaveNumbInvM = 2*Pi*(1.e+06)/LambdaMic

Variable nx = DimSize($InE, 1), nz = DimSize($InE, 2)

Variable xRange = DimDelta($InE, 1)*(nx - 1), zRange = DimDelta($InE, 2)*(nz - 1)
Variable xc = DimOffset($InE, 1) + 0.5*xRange, zc = DimOffset($InE, 2) + 0.5*zRange
Variable xc_mm = 1000.*xc, zc_mm = 1000.*zc, xRange_mm = 1000.*xRange, zRange_mm = 1000.*zRange

//SrwOptThinGenTempl(OutTName,PhEn_keV,xc_mm,zc_mm,xRange_mm,zRange_mm,nx,nz,1)
SrwOptThinTempl(OutTName,1,xc_mm,zc_mm,xRange_mm,zRange_mm,nx,nz)

OutTName += SrwOptThinGenWaveType

$OutTName = cmplx(1, imag(r2polar(srwMiscETransCompens($InE[0](x)(y), WaveNumbInvM, Lsrc, x, y, x0src, z0src, Ph0)))/WaveNumbInvM)
end

//-----------------------------------------------
function/C srwMiscETransCompens(InE, WaveNumb, L, x, z, x0, z0, Ph0)
Variable/C InE
Variable WaveNumb, L, x, z, x0, z0, Ph0

Variable xr = x - x0, zr = z - z0
//Variable Ph = WaveNumb*sqrt(L*L + xr*xr + zr*zr)
Variable Ph = WaveNumb*(xr*xr + zr*zr)/(2*L) + Ph0
//Variable Ph = (xr*xr + zr*zr)/(2*L) + Ph0/WaveNumb // optical path

Variable/C ExpMult = cmplx(cos(Ph), sin(Ph))
if(InE != 0)
	return ExpMult/InE
else
	return ExpMult
endif
end

//-----------------------------------------------
function srwMiscReEAtObs(InE, WaveNumb, L, x, z, x0, z0)
Wave/C InE
Variable WaveNumb, L, x, z, x0, z0

Variable xr = x - x0, zr = z - z0
//Variable Ph = WaveNumb*sqrt(L*L + xr*xr + zr*zr)
Variable Ph = WaveNumb*(xr*xr + zr*zr)/(2*L)
Variable ReMult = cos(Ph), ImMult = sin(Ph)
Variable/C Mult = cmplx(ReMult, ImMult)
Variable/C EAtObs = (InE[0](x)(z))*Mult

return real(EAtObs)
end

//-----------------------------------------------
function srwMiscImEAtObs(InE, WaveNumb, L, x, z, x0, z0)
Wave/C InE
Variable WaveNumb, L, x, z, x0, z0

Variable xr = x - x0, zr = z - z0
//Variable Ph = WaveNumb*sqrt(L*L + xr*xr + zr*zr)
Variable Ph = WaveNumb*(xr*xr + zr*zr)/(2*L)
Variable ReMult = cos(Ph), ImMult = sin(Ph)
Variable/C Mult = cmplx(ReMult, ImMult)
Variable/C EAtObs = (InE[0](x)(z))*Mult

return imag(EAtObs)
end
