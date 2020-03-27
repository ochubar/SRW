#pragma rtGlobals=1		// Use modern global access method.

//==============================================================================
//Calculating Multi-Electron Intensity Distributions and deducing Brightness by Wavefront Propagation
//==============================================================================
proc SrwWfrPropEstimUndBright(nmCoreBright, nmElec, nmMag, nmObsList, nmMagFact, harmStart, harmEnd, nPartInt, nPartProp, sampFact)
string nmCoreBright
string nmElec 
string nmMag
string nmMagFact
string nmObsList //List of 3 names of Observation structures: for on-axis spectra [0], for intensity calc. [1], and for wavefront propag. [2]
variable harmStart
variable harmEnd
variable nPartInt
variable nPartProp
variable sampFact

//Extra parameters (to edit)
variable intMeth = 2 //2- undulator, 3- wiggler
variable intPrec = 0.009
variable intPrecPer = 2
variable relMinPhEn = 0.6
variable relMaxPhEn = 1.4
variable resRange = 1.5
variable resResol = 4
variable resRangeAfterX = 2.5
//End of extra parameters

string nmAuxElec = "auxElecEstimBright_ebm"
string nmAuxElecEnSpOnly = "auxElecEnSpOnlyEstimBright_ebm"
string nmAuxMagCore = "auxMagEstimBright"
string nmAuxMag = nmAuxMagCore + "_mag"
string nmAuxMagBX = nmAuxMagCore + "BX_fld"
string nmAuxMagBZ = nmAuxMagCore + "BZ_fld"
string nmAuxMagPerCore = "auxMagPerEstimBright"
string nmAuxMagPer = nmAuxMagPerCore + "_map"
string nmAuxObsInt = "auxObsIntEstimBright_obs"
string nmAuxObsWfr = "auxObsWfrEstimBright_obs"
string nmAuxSpecCore = "auxWfrSpecEstimBright"
string nmAuxSpecSingleE = nmAuxSpecCore + "_rad"
string nmAuxSpecMultiE = nmAuxSpecCore + "_ras"
string nmAuxSpecIntSuf = "Is"
string nmAuxSpecInt = nmAuxSpecCore + nmAuxSpecIntSuf + "_e"
string nmAuxSpecFluxSuf = "Fm"
string nmAuxSpecFlux = nmAuxSpecCore + nmAuxSpecFluxSuf + "_e"
string nmAuxElecDriftCore = "auxDriftElecEstimBright"
string nmAuxElecDrift = nmAuxElecDriftCore + "_bli"
string nmAuxRadDriftCore = "auxDriftRadEstimBright"
string nmAuxRadDrift = nmAuxRadDriftCore + "_bli"
string nmAuxOptElems = "auxOptElemsList"

string nmResBrightPhotEnCore = nmCoreBright + "EB"
string nmResFluxPhotEnCore = nmCoreBright + "EF"
string nmResBrightCore = nmCoreBright + "B"
string nmResFluxBrightCore = nmCoreBright + "FB"
string nmResFluxCore = nmCoreBright + "F"
string nmResSourceSizeCore = nmCoreBright + "S"
string nmResSourceDivCore = nmCoreBright + "D"
string nmResIntMultiElecCore = nmCoreBright + "I"
string nmResIntWaistMultiElecCore = nmCoreBright + "IW"
string nmResCurBrightPhotEn, nmResCurFluxPhotEn, nmResCurBright, nmResCurFluxBright, nmResCurFlux
string nmResCurSourceSizeX, nmResCurSourceSizeY, nmResCurSourceDivX, nmResCurSourceDivY
string nmResCurIntMultiElecCore, nmResCurIntMultiElec, nmResCurIntWaistMultiElecCore, nmResCurIntWaistMultiElec, nmResCurIntWaistMultiElecInf
string nmResCurIntMultiElecX, nmResCurIntMultiElecY, nmResCurIntWaistMultiElecX, nmResCurIntWaistMultiElecY

duplicate/O $nmElec $nmAuxElec
duplicate/O $nmElec $nmAuxElecEnSpOnly
duplicate/O $($nmObsList[1]) $nmAuxObsInt
duplicate/O $($nmObsList[2]) $nmAuxObsWfr
SrwMagDupl(nmMag, nmAuxMagCore)

variable sStart = dimoffset($nmAuxMagBZ, 0)
variable sStep = dimdelta($nmAuxMagBZ, 0)
variable ns = dimsize($nmAuxMagBZ, 0)
variable sRange = sStep*(ns - 1)
variable maxUndPer = 0.25*sRange*1000

variable elecEn = srwGetElecBeamEnergy(nmElec)
variable elecS0 = srwGetElecBeamLongPos(nmElec)
variable elecSigX = srwGetElecBeamHorSizeRMS(nmElec) //[m]
variable elecSigY = srwGetElecBeamVertSizeRMS(nmElec) //[m]

variable elecDriftLen = (sStart + sStep) - elecS0
SrwOptDrift(nmAuxElecDriftCore, elecDriftLen)
srElecBeamPropag($nmAuxElec, $nmAuxElecDrift)
srElecBeamPropag($nmAuxElecEnSpOnly, $nmAuxElecDrift)
//srwSetElecBeamEmitX(nmAuxElecEnSpOnly, 1e-12)
//srwSetElecBeamEmitZ(nmAuxElecEnSpOnly, 1e-12)
//srwSetElecBeamBetaX(nmAuxElecEnSpOnly, 1)
//srwSetElecBeamBetaZ(nmAuxElecEnSpOnly, 1)
//srwSetElecBeamAlphaX(nmAuxElecEnSpOnly, 0)
//srwSetElecBeamAlphaZ(nmAuxElecEnSpOnly, 0)
SrwElecThick(nmAuxElecEnSpOnly, srwGetElecBeamRelEnSprRMS(nmElec), 1e-12, 1e-12, 1, 1, 0, 0, 0, 0)

variable eStartSpec = srwGetSmpPhotEnStart($nmObsList[0])
variable eEndSpec = srwGetSmpPhotEnEnd($nmObsList[0])

variable xcIntMultiElec = 0.5*(srwGetSmpHorPosStart(nmAuxObsInt) + srwGetSmpHorPosEnd(nmAuxObsInt)) //[m]
variable ycIntMultiElec = 0.5*(srwGetSmpVertPosStart(nmAuxObsInt) + srwGetSmpVertPosEnd(nmAuxObsInt)) //[m]
variable rObsIntMultiElec = srwGetSmpLongPos(nmAuxObsInt) //[m]
SrwOptDrift(nmAuxRadDriftCore, -rObsIntMultiElec)

make/O/T/N=(2, 13) $nmAuxOptElems
$nmAuxOptElems[0][0] = nmAuxRadDrift
$nmAuxOptElems[0][1] = "2"
$nmAuxOptElems[0][2] = "2"
$nmAuxOptElems[0][3] = "1"
$nmAuxOptElems[0][4] = "1"
$nmAuxOptElems[0][5] = "1"
$nmAuxOptElems[0][6] = num2str(resRange)
$nmAuxOptElems[0][7] = num2str(resResol)
$nmAuxOptElems[0][8] = num2str(resRange)
$nmAuxOptElems[0][9] = num2str(resResol)
$nmAuxOptElems[0][10] = "0"
$nmAuxOptElems[0][11] = "0"
$nmAuxOptElems[0][12] = "0"
$nmAuxOptElems[1][5] = "1"
$nmAuxOptElems[1][6] = num2str(resRangeAfterX)
$nmAuxOptElems[1][7] = "1"
$nmAuxOptElems[1][8] = "1"
$nmAuxOptElems[1][9] = "1"
$nmAuxOptElems[1][10] = "0"
$nmAuxOptElems[1][11] = "0"
$nmAuxOptElems[1][12] = "0"

variable nBrightPt = dimsize($nmMagFact, 0)
//variable nHarm = harmEnd - harmStart + 1

variable iBrightPt = 0, iHarm, iMagPerHarm
variable factB, effUndPer, nMagPerHarm, KxE2, KyE2, Keff, Beff, harmResonPhotEn, harmPhotEnMaxFlux
variable infFundPhotEn, supFundPhotEn, d_infFundPhotEn, d_supFundPhotEn, approxResonPhotEn
variable nxIntMultiElec, xStartIntMultiElec, xStepIntMultiElec
variable nyIntMultiElec, yStartIntMultiElec, yStepIntMultiElec
variable xcIntMultiElecWaist, ycIntMultiElecWaist, vCurBrightMult
variable multBright = 1.e-12/(4*Pi*Pi)
string nmAuxMagPerHarm, strHarm, strSufHarmPhotEn, strBrightWinName, strFluxWinName, str2exe
make/D/N=4/O W_coef

do
	factB = $nmMagFact[iBrightPt]
	SrwMagDupl(nmMag, nmAuxMagCore)
	$nmAuxMagBX *= factB
	$nmAuxMagBZ *= factB
	
	SrwMagPrec(nmAuxMag, intMeth, intPrec, intPrec, 10000, 1, 0, 0)
	SrwWfrCreate(nmAuxSpecCore, nmAuxElec, nmAuxMag, $nmObsList[0], 1, 1)
	SrwWfr2Int(nmAuxSpecSingleE, nmAuxSpecIntSuf, 7, 1, 1, 1, 1, 0, 0, 1) //on-axis total single-e intensity vs photon energy

	SrwMagArb2Per(nmAuxMagPerCore, nmAuxMag, 0.03, 5, maxUndPer)
	SrwPerStoCreate(nmAuxSpecCore, nmElec, nmAuxMagPer, $nmObsList[0], 1, harmEnd+4, intPrecPer, intPrecPer, 1)
	SrwSto2IntF(nmAuxSpecMultiE, nmAuxSpecFluxSuf, 7, 1, 1, 1, 0, 0, 1) //multi-e spectral flux
	
	effUndPer = str2num($nmAuxMagPer[0]) //undulator period in [m]
	nMagPerHarm = str2num($nmAuxMagPer[5]) 
	KxE2 = 0
	KyE2 = 0
	iMagPerHarm = 0
	do
		nmAuxMagPerHarm = $nmAuxMagPer[6 + iMagPerHarm]
		if($nmAuxMagPerHarm[1] == 1)
			KyE2 += ($nmAuxMagPerHarm[2]/$nmAuxMagPerHarm[0])^2
		else
			KxE2 += ($nmAuxMagPerHarm[2]/$nmAuxMagPerHarm[0])^2	
		endif
		iMagPerHarm += 1
	while(iMagPerHarm < nMagPerHarm)
	Keff = sqrt(KxE2 + KyE2)
	Beff = Keff/(0.0933729e+03*effUndPer)
	approxResonPhotEn = srUtiUndFundPhotEn(Beff, effUndPer, elecEn, 2) // [eV]
	infFundPhotEn = relMinPhEn*approxResonPhotEn
	supFundPhotEn = relMaxPhEn*approxResonPhotEn	
	d_infFundPhotEn  = approxResonPhotEn - infFundPhotEn
	d_supFundPhotEn = supFundPhotEn - approxResonPhotEn
	
	iHarm = harmStart
	do
		strHarm = num2str(iHarm)
		nmResCurBrightPhotEn = nmResBrightPhotEnCore + strHarm
		nmResCurFluxPhotEn = nmResFluxPhotEnCore + strHarm
		nmResCurBright = nmResBrightCore + strHarm
		nmResCurFluxBright = nmResFluxBrightCore + strHarm
		nmResCurFlux = nmResFluxCore + strHarm
		nmResCurSourceSizeX = nmResSourceSizeCore + "X" + strHarm
		nmResCurSourceSizeY = nmResSourceSizeCore + "Y" + strHarm
		nmResCurSourceDivX = nmResSourceDivCore + "X" + strHarm
		nmResCurSourceDivY = nmResSourceDivCore + "Y" + strHarm

		if(iBrightPt == 0)
			make/O/D/N=(nBrightPt) $nmResCurBrightPhotEn, $nmResCurFluxPhotEn, $nmResCurBright, $nmResCurFluxBright, $nmResCurFlux, $nmResCurSourceSizeX, $nmResCurSourceSizeY, $nmResCurSourceDivX, $nmResCurSourceDivY
		endif

		infFundPhotEn = iHarm*approxResonPhotEn - d_infFundPhotEn
		supFundPhotEn = iHarm*approxResonPhotEn + d_supFundPhotEn
		
		if(!((supFundPhotEn < eStartSpec) %| (infFundPhotEn > eEndSpec)))
		
			WaveStats/Q/R=(infFundPhotEn, supFundPhotEn) $nmAuxSpecInt
			if(V_maxloc <= infFundPhotEn)
				print "WARNING: Lower limit reached when searching harmonic No:", iHarm, "(filament e-beam)"
			endif
			if(V_maxloc >= supFundPhotEn)
				print "WARNING: Upper limit reached when searching harmonic No:", iHarm, "(filament e-beam)"
			endif
			harmResonPhotEn = V_maxloc
			//set this value in Obs struct to calculate Intensity and for Wavefront Prop:
			$nmAuxObsInt[5] = harmResonPhotEn
			$nmAuxObsInt[6] = harmResonPhotEn
			$nmAuxObsWfr[5] = harmResonPhotEn
			$nmAuxObsWfr[6] = harmResonPhotEn
			$nmResCurBrightPhotEn[iBrightPt] = harmResonPhotEn
			
			WaveStats/Q/R=(infFundPhotEn, supFundPhotEn) $nmAuxSpecFlux
			if(V_maxloc <= infFundPhotEn)
				print "WARNING: Lower limit reached when searching harmonic No:", iHarm, "(finite-emittance e-beam)"
			endif
			if(V_maxloc >= supFundPhotEn)
				print "WARNING: Upper limit reached when searching harmonic No:", iHarm, "(finite-emittance e-beam)"
			endif
			harmPhotEnMaxFlux = V_maxloc
			$nmResCurFluxPhotEn[iBrightPt] = harmPhotEnMaxFlux
			$nmResCurFluxBright[iBrightPt] = $nmAuxSpecFlux(harmResonPhotEn)
			$nmResCurFlux[iBrightPt] =  V_max //$nmAuxSpecFlux(harmPhotEnMaxFlux)

			//Calculating Multi-Electron Intensity Distribution to estimate Radiation Beam Angular Divergences
			strSufHarmPhotEn = "H" + num2str(iHarm) + "E" + num2str(harmResonPhotEn*0.001)
			nmResCurIntMultiElecCore = nmResIntMultiElecCore + strSufHarmPhotEn
			nmResCurIntMultiElec = nmResCurIntMultiElecCore + "_xz"
			nmResCurIntMultiElecX = nmResCurIntMultiElecCore + "_x"
			nmResCurIntMultiElecY = nmResCurIntMultiElecCore + "_z"
			SrwIntArbMCCreate(nmResCurIntMultiElecCore, nmAuxElec, nmAuxMag, nmAuxObsInt, 7, nPartInt, intMeth, intPrec, 2)
		
			nxIntMultiElec = dimsize($nmResCurIntMultiElec, 0); nyIntMultiElec = dimsize($nmResCurIntMultiElec, 1)
			xStartIntMultiElec = dimoffset($nmResCurIntMultiElec, 0); yStartIntMultiElec = dimoffset($nmResCurIntMultiElec, 1)
			xStepIntMultiElec = dimdelta($nmResCurIntMultiElec, 0); yStepIntMultiElec = dimdelta($nmResCurIntMultiElec, 1)
		
			make/O/N=(nxIntMultiElec) $nmResCurIntMultiElecX
			SetScale/P x xStartIntMultiElec, xStepIntMultiElec, "m", $nmResCurIntMultiElecX
			$nmResCurIntMultiElecX = $nmResCurIntMultiElec(x)(ycIntMultiElec)
		
			make/O/N=(nyIntMultiElec) $nmResCurIntMultiElecY
			SetScale/P x yStartIntMultiElec, yStepIntMultiElec, "m", $nmResCurIntMultiElecY
			$nmResCurIntMultiElecY = $nmResCurIntMultiElec(xcIntMultiElec)(x)
		
			K0 = 0; CurveFit/Q/W=0/H="1000"/NTHR=0 gauss $nmResCurIntMultiElecX
			$nmResCurSourceDivX[iBrightPt] = W_coef[3]/sqrt(2)/rObsIntMultiElec //Horizontal RMS Angular Divergence [rad]
			K0 = 0; CurveFit/Q/W=0/H="1000"/NTHR=0 gauss $nmResCurIntMultiElecY
			$nmResCurSourceDivY[iBrightPt] = W_coef[3]/sqrt(2)/rObsIntMultiElec //Vertical RMS Angular Divergence

			//Calculating Multi-Electron Intensity Distribution at Waist to estimate Radiation Beam Sizes
			nmResCurIntWaistMultiElecCore = nmResIntWaistMultiElecCore + strSufHarmPhotEn //+ "Ires"
			nmResCurIntWaistMultiElec = nmResCurIntWaistMultiElecCore + "Ires_xz"
			nmResCurIntWaistMultiElecInf = nmResCurIntWaistMultiElec + "_inf"
			nmResCurIntWaistMultiElecX = nmResCurIntWaistMultiElecCore + "Ires_x"
			nmResCurIntWaistMultiElecY = nmResCurIntWaistMultiElecCore + "Ires_z"
			SrwWfrEmitPropStokesMultiE(nmResCurIntWaistMultiElecCore, nmAuxElecEnSpOnly, nmAuxMag, nmAuxObsWfr, sampFact, nmAuxOptElems, nPartProp, 8, 0, 0)

			SrwUtiConvWaveWithGaus2D(nmResCurIntWaistMultiElec, elecSigX, elecSigY)
			srUtiSpotInfo($nmResCurIntWaistMultiElec)
			xcIntMultiElecWaist = $nmResCurIntWaistMultiElecInf[1]
			ycIntMultiElecWaist = $nmResCurIntWaistMultiElecInf[2]
			$nmResCurIntWaistMultiElecX = $nmResCurIntWaistMultiElec(x)(ycIntMultiElecWaist)
			$nmResCurIntWaistMultiElecY = $nmResCurIntWaistMultiElec(xcIntMultiElecWaist)(x)

			K0 = 0; CurveFit/Q/W=0/H="1000"/NTHR=0 gauss $nmResCurIntWaistMultiElecX
			$nmResCurSourceSizeX[iBrightPt] = W_coef[3]/sqrt(2) //Horizontal RMS Size [m]
			K0 = 0; CurveFit/Q/W=0/H="1000"/NTHR=0 gauss $nmResCurIntWaistMultiElecY
			$nmResCurSourceSizeY[iBrightPt] = W_coef[3]/sqrt(2) //Vertical RMS Size
			
			vCurBrightMult = multBright/($nmResCurSourceSizeX[iBrightPt]*$nmResCurSourceDivX[iBrightPt]*$nmResCurSourceSizeY[iBrightPt]*$nmResCurSourceDivY[iBrightPt])
			$nmResCurBright[iBrightPt] = vCurBrightMult*$nmResCurFluxBright[iBrightPt]

			if(iBrightPt == 0)
				if(iHarm == harmStart)
					display $nmResCurBright vs $nmResCurBrightPhotEn
					SrwUtiGraphAddFrameAndGrid()
					ModifyGraph log(left)=1; ModifyGraph log=1
					Label bottom "\\Z12Photon Energy"
					Label left "\\Z12Brightness [Ph/s/.1%bw/mm\\S2\\M\\Z12/mr\\S2\\M\\Z12]"
					SrwUtiGraphWindResize(10,500,350,200,0,0)
					strBrightWinName = WinName(0,1)
				
					display $nmResCurFlux vs $nmResCurFluxPhotEn
					SrwUtiGraphAddFrameAndGrid()
					ModifyGraph log(left)=1; ModifyGraph log=1
					Label bottom "\\Z12Photon Energy"
					Label left "\\Z12Flux [Ph/s/.1%bw]"
					SrwUtiGraphWindResize(310,500,350,200,0,0)
					strFluxWinName = WinName(0,1)
				else
					str2exe = "AppendToGraph/W=" + strBrightWinName + " " + nmResCurBright + " vs " + nmResCurBrightPhotEn
					str2exe += ";AppendToGraph/W=" + strFluxWinName + " " + nmResCurFlux + " vs " + nmResCurFluxPhotEn
					execute str2exe
				endif
			endif

			approxResonPhotEn = harmResonPhotEn/iHarm
		else
			print "WARNING: Search fo harmonic No:", iHarm, "was not made, because it is supposedly outside the calculated spectral range"
		endif
		iHarm += 2
	while(iHarm <= harmEnd)
	iBrightPt += 1
while(iBrightPt < nBrightPt)
end

//==============================================================================
//Auxiliary function for Undulator Radiation Brightness estimation
//==============================================================================
function srwBrilAuxSinc(tx, ty, p, powType)
variable tx, ty, p, powType
variable arg = pi*(p + tx*tx + ty*ty)
if(arg == 0)
	return 1
endif
variable res = sin(arg)/arg
if(powType == 2)
	res *= res
endif
return res
end

//==============================================================================
//Calculate Auxiliary Universal function for radial single-electron intensity distribution at Waist or in Far Field
//==============================================================================
proc SrwBrilUndUnivFuncSingleE(nmUnivFunc, waistOrFarField, Rmax, nR, Pmin, Pmax, nP)
string nmUnivFunc = srwUtiGetValS("nmUnivFunc", "wUndHarmUnivFuncWaist", "SrwBrilUndUnivFuncWaistSingleE")
variable waistOrFarField = srwUtiGetValN("waistOrFarField", 1, "SrwBrilUndUnivFuncWaistSingleE")
variable Rmax = srwUtiGetValN("Rmax", 15, "SrwBrilUndUnivFuncWaistSingleE")
variable nR = srwUtiGetValN("nR", 1000, "SrwBrilUndUnivFuncWaistSingleE")
//variable Prange = srwUtiGetValN("Prange", 10, "SrwBrilUndUnivFuncWaistSingleE")
variable Pmin = srwUtiGetValN("Pmin", -5, "SrwBrilUndUnivFuncWaistSingleE")
variable Pmax = srwUtiGetValN("Pmax", 5, "SrwBrilUndUnivFuncWaistSingleE")
variable nP = srwUtiGetValN("nP", 200, "SrwBrilUndUnivFuncWaistSingleE")
prompt nmUnivFunc, "Name for the Single-E UR Univ. Func. data wave"
prompt waistOrFarField, "Type of the Intensity Univ. Func.", popup "At Waist;In Far Field"
prompt Rmax, "Maximal Normalized Radial Coordinate"
prompt nR, "Number of Points over Radial Coordinate"
//prompt Prange, "Range of Rel. Normalized Energy"
prompt Pmin, "Minimal Rel. Normalized Energy"
prompt Pmax, "Maximal Rel. Normalized Energy"
prompt nP, "Number of Points vs Energy"
silent 1         |       Generating data ...
PauseUpdate

variable lenNmUnivFunc = strlen(nmUnivFunc)
if(lenNmUnivFunc <= 0)
	abort "No new wave name has been provided"
else
	if(lenNmUnivFunc > 31)
		abort "New wave name is too long"
	endif
endif

srwUtiSetValS("nmUnivFunc", nmUnivFunc, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("waistOrFarField", waistOrFarField, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("Rmax", Rmax, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("nR", nR, "SrwBrilUndUnivFuncWaistSingleE")
//srwUtiSetValN("Prange", Prange, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("Pmin", Pmin, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("Pmax", Pmax, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("nP", nP, "SrwBrilUndUnivFuncWaistSingleE")

string nmAuxElecFldHarm = "wAuxElecFldHarmUnivFunc"
string nmAuxElecFldHarmFT = nmAuxElecFldHarm + "FT"
//string nmAuxElecFldHarmFTR = nmAuxElecFldHarmFT + "R"

variable twoNpR = 2*nR
make/C/O/N=(twoNpR,twoNpR) $nmAuxElecFldHarm

variable approxRstep = Rmax/nR
variable angStart = -0.5/approxRstep
variable angStep = abs(angStart)/nR
SetScale/P x angStart,angStep,"", $nmAuxElecFldHarm
SetScale/P y angStart,angStep,"", $nmAuxElecFldHarm

make/O/N=(nR, nP) $nmUnivFunc
//variable pStep = Prange/(nP - 1), pStart = -0.5*Prange
variable pStep = (Pmax - Pmin)/(nP - 1) //, pStart = -0.5*Prange

//SetScale/P y pStart,pStep,"", $nmUnivFunc
SetScale/P y Pmin,pStep,"", $nmUnivFunc

variable ip = 0, pCur = Pmin, rStep
do
	if(waistOrFarField == 1) //at Waist
		$nmAuxElecFldHarm = cmplx(srwBrilAuxSinc(x, y, pCur, 1), 0)
		duplicate/O $nmAuxElecFldHarm $nmAuxElecFldHarmFT

		srFFT2D($nmAuxElecFldHarmFT, -1)

		if(ip == 0)
			//SrwUtiWaveDuplTypeChange(nmAuxElecFldHarmFTR, nmAuxElecFldHarmFT, 3)
			rStep = dimdelta($nmAuxElecFldHarmFT, 0)
			SetScale/P x 0,rStep,"", $nmUnivFunc
		endif
		$nmUnivFunc[][ip] = magsqr($nmAuxElecFldHarmFT[nR + p][nR])
	else //in Far Field
		if(ip == 0)
			rStep = 0
			if(nR > 1)
				rStep = Rmax/(nR - 1)
			endif
			SetScale/P x 0,rStep,"", $nmUnivFunc
		endif

		$nmUnivFunc()[ip] = srwBrilAuxSinc(x, 0, pCur, 2)
	endif

	pCur += pStep
	ip += 1
while(ip < nP)
killwaves/Z $nmAuxElecFldHarm, $nmAuxElecFldHarmFT
end

//==============================================================================
//Calculate Auxiliary Universal function for Radial Intensity distribution at Waist, for one value of Energy Spread
//Takes 2D wave; returns 2D wave.
//==============================================================================
proc SrwBrilUndExtrIntWIthEnSpr(nmRadInt, nmUnivFunc, sigT)
string nmRadInt = srwUtiGetValS("nmRadInt", "wIntWithEnSpr", "SrwBrilUndExtrIntWIthEnSpr")
string nmUnivFunc = srwUtiGetValS("nmUnivFunc", "wUndHarmUnivFuncWaist", "SrwBrilUndUnivFuncWaistSingleE")
variable sigT = srwUtiGetValN("sigT", 0.5, "SrwBrilUndExtrIntWIthEnSpr")
prompt nmRadInt, "Name for Radial UR Intensity Distribution"
prompt nmUnivFunc, "Single-E UR \"Universal Function\" data wave", popup "_none_;" + Wavelist("*",";","TEXT:0,DIMS:2")
prompt sigT, "Normalized RMS Relative Electron Energy Spread"
silent 1         |       Calculating Intensity ...
PauseUpdate

if(!exists(nmUnivFunc))
	abort "Input wave was not found"
endif

variable lenNmRadInt = strlen(nmRadInt)
if(lenNmRadInt <= 0)
	abort "No new wave name has been provided"
else
	if(lenNmRadInt > 31)
		abort "New wave name is too long"
	endif
endif

srwUtiSetValS("nmRadInt", nmRadInt, "SrwBrilUndExtrIntWIthEnSpr")
srwUtiSetValS("nmUnivFunc", nmUnivFunc, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("sigT", sigT, "SrwBrilUndExtrIntWIthEnSpr")

variable nR = dimsize($nmUnivFunc, 0)
variable nP = dimsize($nmUnivFunc, 1), startP = dimoffset($nmUnivFunc, 1), stepP = dimdelta($nmUnivFunc, 1)
variable halfRangeP = 0.5*stepP*(nP - 1)
variable endP = startP + 2*halfRangeP

string nmAuxWaveP = "wAuxInt1dWithEnSpr"
//make/O/N=(nP) $nmAuxWaveP
make/O/N=(nP*2) $nmAuxWaveP

SetScale/P x (startP - halfRangeP),stepP,"", $nmAuxWaveP
duplicate/O $nmUnivFunc $nmRadInt

if(sigT <= 0)
	return //no need to convolute
endif

variable twoSigT = 2*sigT
variable iR=0
do
	$nmAuxWaveP = $nmUnivFunc[iR](x)*srwUtiNonZeroInterval(x, startP, endP)
	SrwUtiConvWaveWithGaus1D(nmAuxWaveP, twoSigT)
	$nmRadInt[iR]() = $nmAuxWaveP(y)

	iR += 1
while(iR < nR)

killwaves/Z $nmAuxWaveP
end

//==============================================================================
//Calculate normalized RMS Size of Radial UR Intensity distribution vs 2 variables:
// - normalized relative "Detuning of Photon Energy" (with respect to on-axis resonant value)
// - normalized relative electron Energy Spread
//==============================================================================
proc SrwBrilUndEstimRMSvsEnergy(nmResRMSvsEn, nmUnivFunc, SigTmax, nSigT, meth, fluxPortion, calcTotFlux)
string nmResRMSvsEn = srwUtiGetValS("nmResRMSvsEn", "wRMSvsEn", "SrwBrilUndEstimRMSvsEnergy")
string nmUnivFunc = srwUtiGetValS("nmUnivFunc", "wUndHarmUnivFuncWaist", "SrwBrilUndUnivFuncWaistSingleE")
variable SigTmax = srwUtiGetValN("SigTmax", 5, "SrwBrilUndEstimRMSvsEnergy")
variable nSigT = srwUtiGetValN("nSigT", 100, "SrwBrilUndEstimRMSvsEnergy")
variable meth = srwUtiGetValN("meth", 5, "SrwBrilUndEstimRMSvsEnergy")
variable fluxPortion = srwUtiGetValN("fluxPortion", 0.9, "SrwBrilUndEstimRMSvsEnergy")
variable calcTotFlux = srwUtiGetValN("calcTotFlux", 1, "SrwBrilUndEstimRMSvsEnergy")
prompt nmResRMSvsEn, "Name for Radial RMS Size wave"
prompt nmUnivFunc, "Single-E UR Univ. Func. data wave", popup "_none_;" + Wavelist("*",";","TEXT:0,DIMS:2")
prompt SigTmax, "Max. Norm. Electron Energy Spread"
prompt nSigT, "Number of Energy Spread values"
prompt meth, "RMS Size Estimation Method", popup "Gaussian Fit;Statistical Analysis"
prompt fluxPortion, "Flux Portion for Stat. Analysis"
prompt calcTotFlux, "Calculate Total Flux as well?", popup "Yes;No"
silent 1         |       Calculating ...
PauseUpdate

if(!exists(nmUnivFunc))
	abort "Input wave was not found"
endif

variable lenNmResRMSvsEn = strlen(nmResRMSvsEn)
if(lenNmResRMSvsEn <= 0)
	abort "No new wave name has been provided"
else
	if(lenNmResRMSvsEn > 30)
		abort "New wave name is too long"
	endif
endif

srwUtiSetValS("nmResRMSvsEn", nmResRMSvsEn, "SrwBrilUndEstimRMSvsEnergy")
srwUtiSetValS("nmUnivFunc", nmUnivFunc, "SrwBrilUndUnivFuncWaistSingleE")
srwUtiSetValN("SigTmax", SigTmax, "SrwBrilUndEstimRMSvsEnergy")
srwUtiSetValN("nSigT", nSigT, "SrwBrilUndEstimRMSvsEnergy")
srwUtiSetValN("meth", meth, "SrwBrilUndEstimRMSvsEnergy")
srwUtiSetValN("fluxPortion", fluxPortion, "SrwBrilUndEstimRMSvsEnergy")
srwUtiSetValN("calcTotFlux", calcTotFlux, "SrwBrilUndEstimRMSvsEnergy")

variable stepSigT = 0
if(nSigT > 1)
	stepSigT = SigTmax/(nSigT - 1)
endif

string nmAuxWaveIntVsRVsDetun = "wAuxIntWithEnSprVsRVsDetun"
string nmAuxWaveR = "wAuxIntWithEnSprEstimRMS"
variable nR = dimsize($nmUnivFunc, 0), startR = dimoffset($nmUnivFunc, 0), stepR = dimdelta($nmUnivFunc, 0)
make/O/N=(nR) $nmAuxWaveR
SetScale/P x startR,stepR,"", $nmAuxWaveR

variable nP = dimsize($nmUnivFunc, 1), startP = dimoffset($nmUnivFunc, 1), stepP = dimdelta($nmUnivFunc, 1)
make/O/N=(nP, nSigT) $nmResRMSvsEn
SetScale/P x startP,stepP,"", $nmResRMSvsEn
SetScale/P y 0,stepSigT,"", $nmResRMSvsEn

string nmResFluxVsEn = nmResRMSvsEn + "F"
if(calcTotFlux == 1)
	make/O/N=(nP, nSigT) $nmResFluxVsEn
	SetScale/P x startP,stepP,"", $nmResFluxVsEn
	SetScale/P y 0,stepSigT,"", $nmResFluxVsEn
endif

variable curSigT = 0, iSigT = 0, iP
variable/C auxRMS_Flux
do
	SrwBrilUndExtrIntWIthEnSpr(nmAuxWaveIntVsRVsDetun, nmUnivFunc, curSigT)
	
	iP = 0
	do
		$nmAuxWaveR = $nmAuxWaveIntVsRVsDetun[p][iP]
		
		if(meth == 1)
			$nmResRMSvsEn[iP][iSigT] = srwUtiRadialDistrGausFitRMS($nmAuxWaveR)
			if(calcTotFlux == 1)
				auxRMS_Flux = srwUtiRadialDistrRMS($nmAuxWaveR, 0, 1)
			endif
		endif
		if(meth == 2)
			auxRMS_Flux = srwUtiRadialDistrRMS($nmAuxWaveR, 2, fluxPortion)
			$nmResRMSvsEn[iP][iSigT] = real(auxRMS_Flux)
		endif
		if(calcTotFlux == 1)
			$nmResFluxVsEn[iP][iSigT] = imag(auxRMS_Flux)
		endif
		
		iP += 1
	while(iP < nP)

	curSigT += stepSigT
	iSigT += 1
while(iSigT < nSigT)

killwaves/Z $nmAuxWaveR, $nmAuxWaveIntVsRVsDetun
end

//==============================================================================
//Auxiliary functions
//==============================================================================
function srwBrilUndBessFact(nMIn, bessArg)
variable nMin, bessArg
variable difBess = BessJ(nMin, bessArg) - BessJ(nMin + 1, bessArg)
return difBess*difBess
end

//==============================================================================
function srwBrilUndBessFactExt(n, K1e2, K2e2)
variable n //harmonic number
variable K1e2, K2e2

//variable dKe2 = K1e2 - K2e2
variable bessArg = 0.25*n*(K1e2 - K2e2)/(1 + 0.5*(K1e2 + K2e2))
variable bessIndMin = (n - 1)/2
variable bessIndMax = bessIndMin + 1

variable bess1 = BessJ(bessIndMin, bessArg)
variable bess2 = BessJ(bessIndMax, bessArg)
variable difBess = bess1 - bess2, sumBess = bess1 + bess2
return K1e2*difBess*difBess + K2e2*sumBess*sumBess
end

//==============================================================================
//function srwBrilUndPhotEnDetunCor(dEperE, relEnSpr, KK, nHarm)
//variable dEperE, relEnSpr, KK, nHarm
function srwBrilUndPhotEnDetunCor(dEperE, relEnSpr, K1e2, K2e2, nHarm)
variable dEperE, relEnSpr, K1e2, K2e2, nHarm

//variable fit_y0 = 0.00024858
variable fit_width = 0.63276

//variable auxMult = nHarm*nHarm*KK*KK/(1 + KK*KK/2)/(fit_width*fit_width)
//variable auxMult = nHarm*nHarm*abs(K1e2 - K2e2)/(1 + (K1e2 + K2e2)/2)/(fit_width*fit_width)
variable auxMult = nHarm*nHarm*(K1e2 + K2e2)/(1 + (K1e2 + K2e2)/2)/(fit_width*fit_width)

variable a_sig = auxMult*2*relEnSpr
variable a_sigE2d2 = a_sig*a_sig/2
variable genFact = 0.5 + 0.5*exp(a_sigE2d2)*(1 - erf(sqrt(a_sigE2d2)))
if(dEperE >= 0)
	//return 1
	return genFact
	//OCTEST
	//print "K=", sqrt(K1e2+K2e2), "srwBrilUndPhotEnDetunCor return:", genFact
endif

//variable corArg = nHarm*KK*sqrt(abs(dEperE)/(1 + KK*KK/2))
//variable relArg = corArg/fit_width
//variable res = fit_y0 + exp(-relArg*relArg)

variable relArg = auxMult*dEperE
variable res = exp(relArg)*genFact

//OCTEST
//print "K=", sqrt(K1e2+K2e2), "srwBrilUndPhotEnDetunCor return:", res

return res
end

//==============================================================================
//Load Universal Data waves required for quick brightness calculation
//==============================================================================
proc SrwBrilUndLoadUnivWaves()
//PathInfo SrwUtiPath
//string strPathUnivWavesFolder = S_path + "SRW Procedures:Development:"
string strPathUnivWavesFolder = srwUtiPathStr() + "SRW Procedures:Development:"
string strPathUnivWaveFlux = strPathUnivWavesFolder + "gwSrwBrilUndHarmUnivFlux.ibw"
string strPathUnivWaveDiv = strPathUnivWavesFolder + "gwSrwBrilUndHarmUnivDiv.ibw"
string strPathUnivWaveSize = strPathUnivWavesFolder + "gwSrwBrilUndHarmUnivSize.ibw"
LoadWave/O/Q strPathUnivWaveFlux
LoadWave/O/Q strPathUnivWaveDiv
LoadWave/O/Q strPathUnivWaveSize
end

//==============================================================================
//Estimate UR Harmonic Brightness (Brilliance)
//Test version, taking into account "detuning" over photon energy
//==============================================================================
proc SrwBrilUndHarmEnDet(BrilName,ElecName,MagName,Kmin,Harm,enDetunPar,NbEnpts,Type,Plot)
string BrilName=SrwElecName+SrwUndName
string ElecName=SrwElecName+SrwElecType
string MagName=SrwUndName+SrwUndType
variable Kmin=SrwKmin
variable Harm=SrwBrilHarm
variable enDetunPar=srwUtiGetValN("enDetunPar", 0, "SrwBrilUndHarmDetun")
variable NbEnpts=SrwNbEnpts
variable Type=SrwType
variable Plot=SrwPlot
prompt BrilName,SrwPBrilName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt Kmin,SrwPKmin
prompt Harm,"Harmonic Number" //SrwPBrilHarm
prompt enDetunPar,"Reson. Ph. Energy \"Detuning\" (dE/En)"
prompt NbEnpts,"Number of Photon Energy Values" //SrwPNbEnpts;
//prompt Type,SrwPSrwType,popup "Phot/s/.1%;Phot/s/.1%/mr2;Phot/s/.1%/mr2/mm2"
prompt Type,SrwPSrwType,popup "Flux [Phot/s/.1%];Angular Flux [Phot/s/.1%/mr2];Brilliance [Phot/s/.1%/mr2/mm2];Horizontal RMS Ang. Div.;Vertical RMS Ang. Div.;Horizontal RMS Source Size;Vertical RMS Source Size"
prompt Plot,SrwPPlot,popup "No;Yes"
Silent 1						|	Estimating Brilliance ...
PauseUpdate

Harm=round(Harm*0.500001)*2-1

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwKmin=Kmin
SrwBrilHarm=Harm
srwUtiSetValN("enDetunPar", enDetunPar, "SrwBrilUndHarmDetun")
SrwNbEnpts=NbEnpts
SrwType=Type
SrwPlot=Plot
SrwBrilName=BrilName

string nmUnivFlux = "gwSrwBrilUndHarmUnivFlux"
string nmUnivDiv = "gwSrwBrilUndHarmUnivDiv"
string nmUnivSize = "gwSrwBrilUndHarmUnivSize"
variable univFluxWaveExists = exists(nmUnivFlux)
variable univDivWaveExists = exists(nmUnivDiv)
variable univSizeWaveExists = exists(nmUnivSize)
if((univFluxWaveExists != 1) %| (univDivWaveExists != 1) %| (univSizeWaveExists != 1))
	SrwBrilUndLoadUnivWaves()
	univFluxWaveExists = exists(nmUnivFlux)
	univDivWaveExists = exists(nmUnivDiv)
	univSizeWaveExists = exists(nmUnivSize)
endif

//  Get Kx and Kz 
variable Kz=0,Kx=0,Per,N,Phase=0
variable phX=0, phZ=0

if($($MagName[6])[1]==1)
	Kz = $($MagName[6])[2]
	phX = $($MagName[6])[3]
else
	Kx = $($MagName[6])[2]
	phZ = $($MagName[6])[3]
endif

if(str2num($MagName[5])==2)
	if($($MagName[7])[1]==1)
		Kz = $($MagName[7])[2]
		Phase = Phase+$($MagName[7])[3]
		phX = $($MagName[7])[3]
	else
		Kx = $($MagName[7])[2]
		Phase = Phase-$($MagName[7])[3]
		phZ = $($MagName[7])[3]
	endif
endif
Per=str2num($MagName[0])
N=str2num($MagName[1])/Per
//print "Kx=",Kx,"  Kz=",Kz,"  Period:",Per,"  Number of Periods:",N,"  Phase bw Field Components:",Phase

// Get Electron Beam Parameters
variable en,relEnSpr,cur,sigx,sigpx,sigz,sigpz,mx,mz
en=$ElecName[0]
relEnSpr=$ElecName[13]
cur=$ElecName[1]
sigx=$ElecName[20] //[m^2]
sigpx=$ElecName[22] //[rad^2]
mx=$ElecName[21] //[m]
sigz=$ElecName[23] //[m^2]
sigpz=$ElecName[25] //[rad^2]
mz=$ElecName[24] //[m]

variable mxE2 = mx*mx, mzE2 = mz*mz

string/G Yy,Xx
if(Type==1)
	Yy=SrwFlStr
endif
if(Type==2)
	Yy=SrwAfStr
endif
if(Type==3)
	Yy=SrwBrStr
endif
if(Type==4)
	Yy="DX"
endif
if(Type==5)
	Yy="DY"
endif
if(Type==6)
	Yy="SX"
endif
if(Type==7)
	Yy="SY"
endif

Yy=BrilName+num2str(Harm)+SrwSepStr+Yy
Xx=BrilName+num2str(Harm)+SrwSepStr+SrwEnStr

//variable K2 = Kx*Kx+Kz*Kz
variable d=sqrt(Kx^4+Kz^4+2*Kx^2*Kz^2*sin(Phase))*Harm

variable KxE2 = Kx*Kx, KzE2 = Kz*Kz
variable ph0 = 0.5*atan((KzE2*sin(2*phX) + KxE2*sin(2*phZ))/(KzE2*cos(2*phX) + KxE2*cos(2*phZ)))
variable K1e2 = KzE2*(cos(phX - ph0))^2 + KxE2*(cos(phZ - ph0))^2
variable K2e2 = KzE2*(sin(phX - ph0))^2 + KxE2*(sin(phZ - ph0))^2
variable Ke2 = KxE2 + KzE2
//variable r1e2 = K1e2/Ke2, r2e2 = K2e2/Ke2

variable h1=(Harm-1)/2
//variable h2=(Harm+1)/2
variable L=N*Per,cst
variable convConstSize, convConstDiv, convConstFlux

make/O/N=(NbEnpts)/D/O $Yy $Xx temp
//setscale/I x 1,Kmin^2/K2,"", temp
setscale/I x 1,Kmin^2/Ke2,"", temp
temp=x
//print d,harm,h1,h2,k2

// Photon Energy [eV]
//$Xx=9.5*en*en/Per/(1+temp*K2/2)*Harm
variable constPhotEn = 9.496376
$Xx = (constPhotEn*en*en*Harm/Per/(1+temp*Ke2/2))*(1 + enDetunPar)
SetScale/P y 0,1,"eV", $Xx

variable normDetun = N*Harm*enDetunPar
variable normEnSpr = N*Harm*relEnSpr
variable factDetunAndEnSpr = Pi/2
if(univFluxWaveExists == 1)
	factDetunAndEnSpr = srwUtiInterp2DBilin(normDetun, normEnSpr, $nmUnivFlux)
	//OCTEST
	//print "Harmonic #", Harm, "  factDetunAndEnSpr=", factDetunAndEnSpr
endif

variable invSqrt2 = 1/sqrt(2)
variable factAngDivDetunAndEnSpr = 0.521402 //to check/edit
if(univDivWaveExists == 1)
	factAngDivDetunAndEnSpr = srwUtiInterp2DBilin(normDetun, normEnSpr, $nmUnivDiv)*invSqrt2 //conv. from Radial to Rectangular representation
endif
variable factAngDivDetunAndEnSprE2 = factAngDivDetunAndEnSpr*factAngDivDetunAndEnSpr

variable factSizeDetunAndEnSpr = 0.408113 //to check/edit
if(univSizeWaveExists == 1)
	factSizeDetunAndEnSpr = srwUtiInterp2DBilin(normDetun, normEnSpr, $nmUnivSize)*invSqrt2 //conv. from Radial to Rectangular representation
endif
variable factSizeDetunAndEnSprE2 = factSizeDetunAndEnSpr*factSizeDetunAndEnSpr

// Flux
if((Type==1) %| (Type==2) %| (Type==3))
	//$Yy=1.431E+14*N*cur*temp*k2/(1+temp*k2/2)*Harm*(bessJ(h1,temp*d/(4+temp*2*K2))-bessJ(h2,temp*d/(4+2*temp*K2)))^2
	
	convConstFlux = 4.5546497e+13
	//$Yy = convConstFlux*N*cur*temp*(K2/(1+temp*K2/2))*Harm*srwBrilUndBessFact(h1, temp*d/(4+temp*2*K2))*factDetunAndEnSpr*srwBrilUndPhotEnDetunCor(enDetunPar, relEnSpr, sqrt(temp*K2), Harm)
	$Yy = convConstFlux*N*cur*(Harm/(1+temp*Ke2/2))*srwBrilUndBessFactExt(Harm, temp*K1e2, temp*K2e2)*factDetunAndEnSpr*srwBrilUndPhotEnDetunCor(enDetunPar, relEnSpr, temp*K1e2, temp*K2e2, Harm)
endif

// Angular Flux
if(Type==2)
	//$Yy*=1.744E*N*en^2/1.431*N*harm/(1+temp*k2/2)
	//$Yy/=sqrt(1+sigpx/(12.4/L/$Xx*1e-7))*sqrt(1+sigpz/(12.4/L/$Xx*1e-7))

	convConstDiv = 2*1.239842e-06/L
	//$Yy = (1.744E+14)*N*N*en*en*cur*temp*K2*(Harm/(1+temp*K2/2))^2*srwBrilUndBessFact(h1, temp*d/(4+temp*2*K2))*srwBrilUndPhotEnDetunCor(enDetunPar, relEnSpr, sqrt(temp*K2), Harm)
	$Yy /= (2e+06*Pi)*sqrt((sigpx + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2)*(sigpz + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2))
endif

// Brilliance
if(Type==3)
	cst=(2*pi)^2*1e12
	convConstSize = 0.5*1.239842e-06*L
	convConstDiv = 2*1.239842e-06/L
	
	//$Yy/=cst*sqrt((sigx+12.4/16/pi/pi*L/$Xx*1e-7)*(sigpx+12.4/L/$Xx*1e-7)-mx*mx)
	//$Yy/=sqrt((sigz+12.4/16/pi/pi*L/$Xx*1e-7)*(sigpz+12.4/L/$Xx*1e-7)-mz*mz)
	$Yy /= cst*sqrt((sigx + (convConstSize/$Xx)*factSizeDetunAndEnSprE2)*(sigpx + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2) - mxE2)
	$Yy /= sqrt((sigz + (convConstSize/$Xx)*factSizeDetunAndEnSprE2)*(sigpz + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2) - mzE2)
endif

//Horizontal RMS Ang. Div.
if(Type==4)
	convConstDiv = 2*1.239842e-06/L
	$Yy = sqrt(sigpx + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2)
endif
//Vertical RMS Ang. Div.
if(Type==5)
	convConstDiv = 2*1.239842e-06/L
	$Yy = sqrt(sigpz + (convConstDiv/$Xx)*factAngDivDetunAndEnSprE2)
endif
//Horizontal RMS Source Size
if(Type==6)
	convConstSize = 0.5*1.239842e-06*L
	$Yy = sqrt(sigx + (convConstSize/$Xx)*factSizeDetunAndEnSprE2)
endif
//Vertical RMS Source Size
if(Type==7)
	convConstSize = 0.5*1.239842e-06*L
	$Yy = sqrt(sigz + (convConstSize/$Xx)*factSizeDetunAndEnSprE2)
endif

if(Plot==2)
	Display $Yy vs $Xx
endif
killwaves/Z temp
end 

//==============================================================================
//Estimate UR Brightness (Brilliance)
//Test version, taking into account "detuning" over photon energy
//==============================================================================
proc SrwBrilUndEnDet(BrilName,ElecName,MagName,Kmin,HarmMin,HarmMax,enDetunParH1, NbEnpts,Type,Plot)
string BrilName=SrwElecName+SrwUndName
string ElecName=SrwElecName+SrwElecType
string MagName=SrwUndName+SrwUndType
variable Kmin=SrwKmin
variable HarmMin=SrwBrilHarmMin
variable HarmMax=SrwBrilHarmMax
variable enDetunParH1=srwUtiGetValN("enDetunParH1", 0, "SrwBrilUndEnDet")
variable NbEnpts=SrwNbEnpts
variable Type=SrwType
variable Plot=SrwPlot
prompt BrilName,SrwPBrilName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt Kmin,SrwPKmin
prompt HarmMin,SrwPBrilHarmMin
prompt HarmMax,SrwPBrilHarmMax
prompt enDetunParH1,"Photon Energy \"Detuning\" (dE/E1)"
prompt NbEnpts,SrwPNbEnpts
prompt Type,SrwPSrwType,popup "Flux [Phot/s/.1%];Angular Flux [Phot/s/.1%/mr2];Brilliance [Phot/s/.1%/mr2/mm2];Horizontal RMS Ang. Div.;Vertical RMS Ang. Div.;Horizontal RMS Source Size;Vertical RMS Source Size"
prompt Plot,SrwPPlot,popup "No;Yes"
Silent 1						|	Estimating Brilliance ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if(MagWavePresent == 1)
		SrwBrilUnd()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	SrwBrilUnd()
	Return
endif
if(SrwUndIsEllips(MagName)==0)
	Abort "Sorry, this type of computation supports only sinusoidal magnetic field."
endif
if(HarmMin>HarmMax)
	Abort "Wrong Harmonic Numbers"
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwKmin=Kmin
SrwBrilHarmMin=HarmMin
SrwBrilHarmMax=HarmMax
srwUtiSetValN("enDetunParH1", enDetunParH1, "SrwBrilUndEnDet")
SrwNbEnpts=NbEnpts
SrwType=Type
SrwBrilName=BrilName

string/G Yy,Xx

HarmMin=Max(HarmMin,1)
variable h=HarmMin

if(Plot==2)
	Display 
endif
do
	//SrwBrilUndHarm(BrilName, ElecName, MagName, Kmin,h,NbEnpts,Type,1)
	SrwBrilUndHarmEnDet(BrilName,ElecName,MagName,Kmin,h,enDetunParH1/h,NbEnpts,Type,1)
	
	if(Plot==2)
		AppendToGraph $Yy vs $Xx
	endif
	h+=2
while (h<harmMax+1)

if(Plot==2)
	ModifyGraph log(left)=1
	ModifyGraph log=1
	Label bottom SrwPLabelPhotEn

	if(Type==1)
		Label left SrwPUnitSpAngFlux
	endif
	if(Type==2)
		Label left SrwPUnitSpAngFluxPerUnAngle
	endif
	if(Type==3)
		Label left SrwPUnitBrilliance
	endif
	if((Type==4) %| (Type==5))
		Label left "rad"
	endif
	if((Type==6) %| (Type==7))
		Label left "m"
	endif

	SrwUtiGraphAddFrameAndGrid()
endif

SrwPlot=Plot
end

//===========================================================================
function SrwAuxBrilAddHarm(wF, wE, wFtoAdd, wEtoAdd)
wave wF, wE, wFtoAdd, wEtoAdd
variable np = dimsize(wF, 0), i, en
for(i=0; i<np; i+=1)
	en = wE[i]
	FindLevel/Q wEtoAdd, en
	if(V_flag == 0)
		wF[i] += wFtoAdd[V_LevelX]
	endif
endfor
end

//===========================================================================
function SrwAuxBrilEstimFromFlux(wF, wE, wDX, wDY, wSX, wSY, wDSE)
wave wF, wE, wDX, wDY, wSX, wSY, wDSE
variable cst = (2*pi)^2*1e+12
variable np = dimsize(wF, 0), i, en, i0
for(i=0; i<np; i+=1)
	en = wE[i]
	FindLevel/Q wDSE, en
	if(V_flag == 0)
		i0 = V_LevelX
	else
		if(en <= wDSE[0])
			i0 = 0
		endif
		if(en >= wDSE[dimsize(wDSE, 0) - 1])
			i0 = dimsize(wDSE, 0) - 1
		endif
	endif
	wF[i] /= cst*wDX[i0]*wDY[i0]*wSX[i0]*wSY[i0]
endfor
end

//==============================================================================
//Estimate Apparent Flux / Brightness of an "Adaptive-Gap Undulator"
//Test version, taking into account "detuning" over photon energy
//==============================================================================
//proc SrwBrilUndEnDet(BrilName,ElecName,MagName,Kmin,HarmMin,HarmMax,enDetunParH1, NbEnpts,Type,Plot)
proc SrwBrilAGUEnDet(BrilName,ElecName,nmAguStruct,Kmin,HarmMin,HarmMax,enDetunParRel, NbEnpts,Type,Plot)
string BrilName=srwUtiGetValS("BrilName", "undFlux", "SrwBrilAGUEnDet")
string ElecName=SrwElecName+SrwElecType
string nmAguStruct=srwUtiGetValS("nmAguStruct", "wDescrAGU", "SrwBrilAGUEnDet")
variable Kmin=SrwKmin
variable HarmMin=SrwBrilHarmMin
variable HarmMax=SrwBrilHarmMax
variable enDetunParRel=srwUtiGetValN("enDetunParRel", 0, "SrwBrilAGUEnDet")
variable NbEnpts=SrwNbEnpts
variable Type=SrwType
variable Plot=srwUtiGetValN("Plot", 2, "SrwBrilAGUEnDet")
prompt BrilName,SrwPBrilName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt nmAguStruct,"Name of AGU descr. (per., K, len) wave",popup Wavelist("*",";","TEXT:0,DIMS:2,MINCOLS:3,MAXCOLS:3")
prompt Kmin,SrwPKmin
prompt HarmMin,SrwPBrilHarmMin
prompt HarmMax,SrwPBrilHarmMax
prompt enDetunParRel,"Photon Energy \"Detuning\" (Nu*dE/E1)"
prompt NbEnpts,SrwPNbEnpts
prompt Type,SrwPSrwType,popup "Flux [Phot/s/.1%];Angular Flux [Phot/s/.1%/mr2];Brilliance [Phot/s/.1%/mr2/mm2];Horizontal RMS Ang. Div.;Vertical RMS Ang. Div.;Horizontal RMS Source Size;Vertical RMS Source Size"
prompt Plot,SrwPPlot,popup "No;Yes"
Silent 1						|	Estimating Brilliance ...
PauseUpdate

if(Type > 3)
	abort "Calculation of this radiation characteristic is not supported"
endif

srwUtiSetValS("BrilName", BrilName, "SrwBrilAGUEnDet")
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
//SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
//idbUtiSetValS("nmAguStruct", nmAguStruct, "SrwBrilAGUEnDet")
srwUtiSetValS("nmAguStruct", nmAguStruct, "SrwBrilAGUEnDet") //OC16122019
SrwKmin=Kmin
SrwBrilHarmMin=HarmMin
SrwBrilHarmMax=HarmMax
srwUtiSetValN("enDetunParRel", enDetunParRel, "SrwBrilAGUEnDet")
SrwNbEnpts=NbEnpts
SrwType=Type
SrwBrilName=BrilName
srwUtiSetValN("Plot", Plot, "SrwBrilAGUEnDet")

//variable FluxOrBright = 2 //1- flux, 2- brightness
//string nmAguStruct = "aguBetY2_86SectPerK"
//string nmElec = "NSLSII_HB_fin_tune_ebm"

string nmResFluxCore = BrilName
string nmUndCore = "aguSegmAuxBril"
///variable Kmin = 0.2
//variable nHmax = 31
//variable enDetunPar = -0.4 //-1
//variable nPt = 1000

SrwAllowPrintingExtraInfo = 0

variable nIndepSect = dimsize($nmAguStruct, 0)
//variable sStep = dimdelta($nmAguStruct, 0)
//variable s = dimoffset($nmAguStruct, 0)
variable iSect=0, curPer, curK, curL, iHarm, avgPer=0, avgK=0, totL=0
string nmCurF, nmCurE, nmResF, nmResE, nmUnd, nmDy, nmSy, nmDx, nmSx, nmDSE, nmResB

do
	//if(iSect == 0)
	//	curL = 2*s
	//else
	//	curL = 2*sStep
	//endif
	
	nmUnd = nmUndCore + num2str(iSect)
	curPer = $nmAguStruct[iSect][0]*1000
	curK = $nmAguStruct[iSect][1]
	curL = $nmAguStruct[iSect][2]
	SrwMagPerCreate2D(nmUnd, curPer, curK, 0, curL, 0,1,0,0)

	//calc flux from current segment (without display)
	//SrwBrilUndEnDet(nmUnd, nmElec, nmUnd + "_map", Kmin, 1, nHmax, enDetunParRel/(curL*1000/curPer), nPt, 1, 1)
	SrwBrilUndEnDet(nmUnd, ElecName, nmUnd + "_map", Kmin, HarmMin, HarmMax, enDetunParRel/(curL*1000/curPer), NbEnpts, 1, 1)
	
	iHarm = HarmMin
	do
		nmCurF = nmUnd + num2str(iHarm) + ".F"
		nmCurE = nmUnd +  num2str(iHarm) + ".E"
		
		nmResF = nmResFluxCore + num2str(iHarm) + "F"
		nmResE = nmResFluxCore + num2str(iHarm) + "E"
		
		if(iSect == 0)
			duplicate/O $nmCurF $nmResF
			duplicate/O $nmCurE $nmResE

			if(Plot == 2)
				if(iHarm == HarmMin)
					display $nmResF vs $nmResE
					ModifyGraph log(left)=1
					ModifyGraph log=1
					SrwUtiGraphAddFrameAndGrid()
				else
					appendtograph $nmResF vs $nmResE
				endif
			endif
		else
			SrwAuxBrilAddHarm($nmResF, $nmResE, $nmCurF, $nmCurE)
		endif
	
		iHarm += 2
	while(iHarm <= HarmMax)

	avgPer += curPer
	avgK += curK
	totL += curL

	iSect += 1
while(iSect < nIndepSect)
avgPer /= nIndepSect
avgK /= nIndepSect

if(Type == 3) //Brightness needed

	SrwMagPerCreate2D(nmUndCore + "DS", avgPer, avgK, 0, totL, 0,1,0,0)

	SrwBrilUndEnDet(nmResFluxCore + "DS", ElecName, nmUndCore + "DS_map",Kmin, HarmMin, HarmMax, enDetunParRel/(totL*1000/avgPer), NbEnpts, 4, 1) //Dx
	SrwBrilUndEnDet(nmResFluxCore + "DS", ElecName, nmUndCore + "DS_map",Kmin, HarmMin, HarmMax, enDetunParRel/(totL*1000/avgPer), NbEnpts, 5, 1) //Dy
	SrwBrilUndEnDet(nmResFluxCore + "DS", ElecName, nmUndCore + "DS_map",Kmin, HarmMin, HarmMax, enDetunParRel/(totL*1000/avgPer), NbEnpts, 6, 1) //Sx
	SrwBrilUndEnDet(nmResFluxCore + "DS", ElecName, nmUndCore + "DS_map",Kmin, HarmMin, HarmMax, enDetunParRel/(totL*1000/avgPer), NbEnpts, 7, 1) //Sy

	iHarm = HarmMin
	do
		nmDx = nmResFluxCore + "DS" + num2str(iHarm) + ".DX"
		nmDy = nmResFluxCore + "DS" + num2str(iHarm) + ".DY"
		nmSx = nmResFluxCore + "DS" + num2str(iHarm) + ".SX"
		nmSy = nmResFluxCore + "DS" + num2str(iHarm) + ".SY"
		nmDSE = nmResFluxCore + "DS" + num2str(iHarm) + ".E"
		
		nmResF = nmResFluxCore + num2str(iHarm) + "F"
		nmResB = nmResFluxCore + num2str(iHarm) + "B"
		duplicate/O $nmResF $nmResB
		nmResE = nmResFluxCore + num2str(iHarm) + "E"

		SrwAuxBrilEstimFromFlux($nmResB, $nmResE, $nmDx, $nmDy, $nmSx, $nmSy, $nmDSE)
		
		if(Plot == 2)
			if(iHarm == HarmMin)
				display $nmResB vs $nmResE
				ModifyGraph log(left)=1
				ModifyGraph log=1
				SrwUtiGraphAddFrameAndGrid()
			else
				appendtograph $nmResB vs $nmResE
			endif
		endif
		iHarm += 2
	while(iHarm <= HarmMax)
endif
SrwType=Type
end