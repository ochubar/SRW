
//+++++++++++++++++++++++++++++++++++++++
//
//Function for determining Electron Energy Spread by fitting
//
//+++++++++++++++++++++++++++++++++++++++
function srwElecDiagEnSprFromUndSpecFunc(wParams, enSpr)
wave wParams
variable enSpr

SVAR nmElec_DiagEnSpr, nmUnd_DiagEnSpr, nmObs_DiagEnSpr, nmSpecMeas_DiagEnSpr

srwSetElecBeamRelEnSprRMS(nmElec_DiagEnSpr, enSpr)
string nmCoreStokes = "AuxStokesOptim"
string nmCalcSpecFlux = nmCoreStokes + "F_e"

variable cenHarm = wParams[0]
variable numNeibHarm = wParams[1]
variable precFactStokes = wParams[2]
variable startPhotEn_eV = wParams[3]
variable endPhotEn_eV = wParams[4]

if(cenHarm < 1)
	cenHarm = 1
endif
if(numNeibHarm < 1)
	numNeibHarm = 1
endif
variable startHarm = cenHarm - numNeibHarm, endHarm = cenHarm + numNeibHarm
if(startHarm < 1)
	startHarm = 1
endif
string sPrec = num2str(precFactStokes)
string sToExe = "SrwPerStoCreate(\"" + nmCoreStokes + "\",\"" + nmElec_DiagEnSpr + "\",\"" + nmUnd_DiagEnSpr + "\",\"" + nmObs_DiagEnSpr + "\"," + num2str(startHarm) + "," + num2str(endHarm) + "," + sPrec + "," + sPrec + ",1)"
execute sToExe
sToExe = "SrwSto2IntF(\"" + nmCoreStokes + "_ras\",\"F\",7,1,1,1.,2.5e-10,2.5e-10,1)"
execute sToExe

wave wSpecCalc = $nmCalcSpecFlux
wave wSpecMeas = $nmSpecMeas_DiagEnSpr
wavestats/Q/R=(startPhotEn_eV,endPhotEn_eV) wSpecCalc
variable calcMaxLoc = V_maxloc, calcMaxVal = V_max
wavestats/Q/R=(startPhotEn_eV,endPhotEn_eV) wSpecMeas
variable measMaxLoc = V_maxloc, measMaxVal = V_max

KillVariables/Z V_minloc

wSpecCalc *= (measMaxVal/calcMaxVal)
variable calcExtCoef = measMaxLoc/calcMaxLoc
//variable oldStartArgCalc = dimoffset(wSpecCalc, 0), oldStepArgCalc = dimdelta(wSpecCalc, 0)
//variable indMaxArgCalc = trunc((calcMaxLoc - oldStartArgCalc)/oldStepArgCalc + 0.00001)
//variable newStartArgCalc = oldStartArgCalc*calcExtCoef
//variable newStepArgCalc = oldStepArgCalc
//if(indMaxArgCalc > 0)
//	newStepArgCalc = calcExtCoef*(calcMaxLoc - oldStartArgCalc)/indMaxArgCalc
//endif
//SetScale/P x newStartArgCalc,newStepArgCalc,"", wSpecCalc

//variable leftExtCoef = startPhotEn_eV/measMaxLoc, rightExtCoef = endPhotEn_eV/measMaxLoc
//variable absDifMaxLoc = abs(calcMaxLoc - measMaxLoc)
//variable leftExtCoef = (measMaxLoc - 0.05*absDifMaxLoc)/measMaxLoc, rightExtCoef = (measMaxLoc + 0.05*absDifMaxLoc)/measMaxLoc
variable leftExtCoef = calcExtCoef*(1-0.03), rightExtCoef = calcExtCoef*(1+0.03)

variable tolExtCoef = 0.000001
make/O wParamExtOpt = {startPhotEn_eV, endPhotEn_eV}
string/G nmSpecCalcMod_DiagEnSpr = nmCalcSpecFlux
variable/G gCurWorkMinFitVal = 1e+23, gCurWorkSpecExtCoef = 1
Optimize/A=0/H=(rightExtCoef)/L=(leftExtCoef)/T=(tolExtCoef) srwElecDiagExtCoefFunc, wParamExtOpt
//Optimize/H=(rightExtCoef)/L=(leftExtCoef)/I=(30) srwElecDiagExtCoefFunc, wParamExtOpt

variable/G gCurBestEnSpr, gCurMinFitVal, gCurBestSpecExtCoef
//if(gCurMinFitVal > V_min)
if(gCurMinFitVal > gCurWorkMinFitVal)
	gCurMinFitVal = gCurWorkMinFitVal
	gCurBestSpecExtCoef = gCurWorkSpecExtCoef
	gCurBestEnSpr = enSpr
endif

return V_min
end

//+++++++++++++++++++++++++++++++++++++++
//
//Function for determining spectrum extension coefficient
//
//+++++++++++++++++++++++++++++++++++++++
function srwElecDiagExtCoefFunc(wParams, extCoef)
wave wParams
variable extCoef

SVAR nmSpecMeas_DiagEnSpr, nmSpecCalcMod_DiagEnSpr
string nmAuxFitWave = "wAuxElecDiagEnSprFromUndSpec"

wave wSpecCalc = $nmSpecCalcMod_DiagEnSpr
wave wSpecMeas = $nmSpecMeas_DiagEnSpr

variable oldStartArgCalc = dimoffset(wSpecCalc, 0), oldStepArgCalc = dimdelta(wSpecCalc, 0)
variable npArgCalc = dimsize(wSpecCalc, 0)
variable oldEndArgCalc = oldStartArgCalc + oldStepArgCalc*(npArgCalc - 1)
variable newStartArgCalc = oldStartArgCalc*extCoef
variable newEndArgCalc = oldEndArgCalc*extCoef

duplicate/O $nmSpecCalcMod_DiagEnSpr $nmAuxFitWave
wave wAuxFit = $nmAuxFitWave
SetScale/I x newStartArgCalc,newEndArgCalc,"", wAuxFit

variable startPhotEn_eV = wParams[0]
variable endPhotEn_eV = wParams[1]

//if(exists(nmAuxFitWave) == 0)
//	duplicate/O $nmSpecMeas_DiagEnSpr $nmAuxFitWave
//endif

//wAuxFit = 0
//wAuxFit = srwUtiPow2(wSpecCalc(x) - wSpecMeas(x))*srwUtiNonZeroInterval(x, startPhotEn_eV, endPhotEn_eV)
//wavestats/Q/R=(startPhotEn_eV,endPhotEn_eV) wAuxFit
//variable curMinVal = sqrt(V_Sum)

variable sumDifE2 = 0, x = startPhotEn_eV, xStep = dimdelta(wSpecMeas, 0)
variable b1 = 0
do
	b1 = wAuxFit(x) - wSpecMeas(x)
	sumDifE2 += b1*b1
	x += xStep
while(x < endPhotEn_eV)
variable curMinVal = sqrt(sumDifE2)

	//DoUpdate

KillVariables/Z V_minloc

variable/G gCurWorkMinFitVal, gCurWorkSpecExtCoef
if(gCurWorkMinFitVal > curMinVal)
	gCurWorkMinFitVal = curMinVal
	gCurWorkSpecExtCoef = extCoef
endif

return curMinVal
end

//+++++++++++++++++++++++++++++++++++++++
//
//Proc for determining Electron Energy Spread by fitting - TO IMPROVE!
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecDiagEnSprFromUndSpec() //to implement

string/G nmElec_DiagEnSpr = "SOL_SSopt_ebm"
string/G nmUnd_DiagEnSpr = "MagU20_G9_DANPerTm_map"
string/G nmObs_DiagEnSpr = "ObsEPX1s_obs"
string/G nmSpecMeas_DiagEnSpr = "wSpecExp2_proc_nd"

variable cenHarm = 3
variable numNeibHarm = 2
variable precFactStokes = 5
variable startPhotEn_eV = 7450 //to calculate form und. and e-beam
variable endPhotEn_eV = 7700 //to calculate form und. and e-beam
make/O wParamOptEnSpr = {cenHarm, numNeibHarm, precFactStokes, startPhotEn_eV, endPhotEn_eV}
variable lowerEnSpr = 0.0003
variable upperEnSpr = 0.002
variable tolEnSpr = 0.000001

variable/G gCurMinFitVal = 1e+23, gCurBestEnSpr, gCurBestSpecExtCoef

Optimize/Q/L=(lowerEnSpr)/H=(upperEnSpr)/T=(tolEnSpr) srwElecDiagEnSprFromUndSpecFunc, wParamOptEnSpr
print "Rel. Energy Spread found:", gCurBestEnSpr
print "Spec. Extension Coef.:", gCurBestSpecExtCoef
end

//+++++++++++++++++++++++++++++++++++++++
//
//Function for determining projected electron size (hor.) by fitting
//
//+++++++++++++++++++++++++++++++++++++++
function srwElecDiagProjSizeFromIntFunc(wParams, projSizeRMS, multFact)
wave wParams
variable projSizeRMS, multFact


	//return (projSizeRMS - 3)*(projSizeRMS - 3) + (multFact - 1)*(multFact - 1)

SVAR nmIntFilamCalc_DiagProjSize, nmIntMeas_DiagProjSize

string nmAuxWave = "wAuxElecDiagProjSizeFromIntFunc"
if(exists(nmAuxWave) == 0)
	duplicate/O $nmIntFilamCalc_DiagProjSize $nmAuxWave
endif
wave wAux = $nmAuxWave
wave wIntFilamCalc = $nmIntFilamCalc_DiagProjSize
wave wIntMeas = $nmIntMeas_DiagProjSize
wAux = wIntFilamCalc

variable xStart = wParams[0], xEnd = wParams[1]

string sExe = "SrwUtiConvWaveWithGaus1D(\"" + nmAuxWave + "\"," + num2str(projSizeRMS) + ")"
execute sExe
wavestats/Q/R=(xStart, xEnd) wAux
variable curMaxVal = V_max
wAux *= (multFact/curMaxVal)

	DoUpdate

variable sumDifE2 = 0, x = xStart, xStep = dimdelta(wAux, 0)
variable b1 = 0
do
	b1 = wAux(x) - wIntMeas(x)
	sumDifE2 += b1*b1
	x += xStep
while(x < xEnd)
variable np = trunc((xEnd - xStart)/xStep + 0.00001) + 1
variable curMinVal = sqrt(sumDifE2/np)

variable/G gCurMinFitVal, gCurBestProjSize, gCurBestMultFact
if(gCurMinFitVal > curMinVal)
	gCurMinFitVal = curMinVal
	gCurBestProjSize = projSizeRMS
	gCurBestMultFact = multFact
endif
return curMinVal
end

//+++++++++++++++++++++++++++++++++++++++
//
//Proc for determining projectted RMS size by fitting - TO IMPROVE!
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecDiagProjSizeFromInt() //to implement

string/G nmIntFilamCalc_DiagProjSize = "wAuxCalcIntCutHs"
string/G nmIntMeas_DiagProjSize = "wMeasCutHsm"

variable xStart = 0 //-5e-03
variable xEnd = 6e-03
make/O wParamOptProjSize = {xStart, xEnd}
//variable lowerProjSizeRMS = 0.0001 //[m]
//variable upperProjSizeRMS = 0.001 //[m]
//variable tolProjSizeRMS = 0.000001
variable expectProjSizeRMS = 0.0006
variable nDigitsCor = 6

wavestats/Q/R=(xStart, xEnd) $nmIntMeas_DiagProjSize
variable expectMultFact = V_max

variable/G gCurMinFitVal = 1e+23, gCurBestProjSize, gCurBestMultFact
//Optimize/Q/D=(nDigitsCor)/R = {expectProjSizeRMS, expectMultFact}/X = {expectProjSizeRMS, expectMultFact}/M = {0, 0} srwElecDiagProjSizeFromIntFunc, wParamOptProjSize
//Optimize/Q/X = {expectProjSizeRMS, expectMultFact} srwElecDiagProjSizeFromIntFunc, wParamOptProjSize

Optimize/Q/D=4/R = {0.0006, 200}/X = {0.0006, 200}/M = {1, 0} srwElecDiagProjSizeFromIntFunc, wParamOptProjSize

//optimize/D=4/M={1,0}/X={5, 10} srwElecDiagEnSprFromUndSpec, wAuxParams

print "Projected RMS Size found:", gCurBestProjSize
print "Corresponding Multiplication Factor:", gCurBestMultFact
end