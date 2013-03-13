//#pragma rtGlobals=1		// Use modern global access method

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates Fundamental Photon Energy and Polarization Rate vs Two Undulator Parameters 
//(e.g. Gap and Phase) from (measured) Magnetic Field.
//Call example: 
//SrwUndPhotEnPolarVsParam("otHU80_PLEIADES", "root:Hall:Map:", "otNameCoreMeas_B_vs_S", "SOL_MSoff_ebm", "ObsE_obs", 1, 0.5, 1.5, 2, 0.005)
//
//First prepare a wave with names of waves containing magnetic measurements data:
//make/T/O/N=(13,17,2) otMeasNamesBvsS
//edit otMeasNamesBvsS
//make/O/N=13 otMeasGapValues
//edit otMeasGapValues
//make/O/N=17 otMeasPhaseValues
//appendtotable otMeasPhaseValues
//otMeasNamesBvsS[][][0] = "Bx_Mes_G" + num2str(otMeasGapValues[p]*10) + "_P" + ReplaceString("-", num2str(otMeasPhaseValues[q]*10), "m") + "_II_C"
//otMeasNamesBvsS[][][1] = "Bz_Mes_G" + num2str(otMeasGapValues[p]*10) + "_P" + ReplaceString("-", num2str(otMeasPhaseValues[q]*10), "m") + "_II_C"

//Inversion of these tables:
//make/O/N=(10,17) otFundPhotEnVsParam, otCircPolRateVsParam
//otFundPhotEnVsParam = otHU80_PLEIADESPhEn[p][q][0]
//otCircPolRateVsParam = otHU80_PLEIADESPolR[p][q][1][0]
//
//FuncFit/Q srwUndPhotEnVsGap W_coef  wAuxColTable /X=otGapValues /W=otGapValuesWeights /I=1 /D
//edit otWcoef01PhotEn
//SrwUtiCalcFitPolyCoef2D("otPhotEnTestCoefs","otFundPhotEnVsParam","otGapValues","otPhaseValues","srwUndPhotEnVsGap", "srwUndPhotEnVsPhase","otGapValuesWeights","otPhaseValuesWeights","otWcoef01PhotEn","otWcoef02PhotEn")
//SrwUtiCalcFitPolyCoef2D("otPolRateTestCoefs","otCircPolRateVsParam","otGapValues","otPhaseValues","srwUndPolRateVsGap","srwUndPolRateVsPhase","otGapValuesWeights","otPhaseValuesWeights","otWcoef01PolRate","otWcoef02PolRate")
//
//duplicate/O otGapVsFundPhotEn_HU80LH otGapValuesDense
//otGapValuesDense = 32 + p*2
//SrwUndFindGapAndPhaseFromParam("otGapVsFundPhotEn_HU80CR", "otPhaseVsFundPhotEn_HU80CR", "otGapValuesDense", 1.1)
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUndPhotEnPolarVsParam(nmCoreRes, nmFolderMeas, nmCoreNamesMeas, nmElec, nmObsList, numHarm, relMinPhEn, relMaxPhEn, intPrec, rangePer)
//string nmResPhotEnVsParam = srwUtiGetValS("nmResPhotEnVsParam","otFundPhotEnVsParam","SrwUndFundPhotEnPolarVsParam")
//string nmResPolRateVsParam = srwUtiGetValS("nmResPolRateVsParam","otPolRateVsParam","SrwUndFundPhotEnPolarVsParam")
string nmCoreRes = srwUtiGetValS("nmCoreRes","otUnd","SrwUndPhotEnPolarVsParam")
//variable approxUndPer = srwUtiGetValN("approxUndPer",80,"SrwUndFundPhotEnPolarVsParam")
string nmFolderMeas = srwUtiGetValS("nmFolderMeas","root:Hall:","SrwUndPhotEnPolarVsParam")
string nmCoreNamesMeas = srwUtiGetValS("nmCoreNamesMeas","otNameCore_B_vs_S","SrwUndPhotEnPolarVsParam")
//string sufNamesMeas = srwUtiGetValS("sufNamesMeas","II","SrwUndPhotEnPolarVsParam")
string nmElec = srwUtiGetValS("SrwElecName","Elec","")+SrwElecType
//string nmObs = srwUtiGetValS("SrwSmpName","Obs","")+SrwSmpType
string nmObsList = srwUtiGetValS("nmObsList","ObsList","SrwUndPhotEnPolarVsParam")
variable numHarm = srwUtiGetValN("numHarm",1,"SrwUndPhotEnPolarVsParam")
variable relMinPhEn = srwUtiGetValN("relMinPhEn",0.95,"SrwUndPhotEnPolarVsParam")
variable relMaxPhEn = srwUtiGetValN("relMaxPhEn",1.5,"SrwUndPhotEnPolarVsParam")
//variable intMeth = srwUtiGetValN("intMeth",3,"SrwUndPhotEnPolarVsParam")
variable intPrec = srwUtiGetValN("intPrec",0.007,"SrwUndPhotEnPolarVsParam")
variable rangePer = srwUtiGetValN("rangePer",0,"SrwUndPhotEnPolarVsParam")
//prompt nmResPhotEnVsParam,"Name for Photon En. data wave to calc."
//prompt nmResPolRateVsParam,"Name for Polar. Rate data wave to calc."
//prompt approxUndPer,"Approximate Undulator Period [mm]"
prompt nmCoreRes,"Core Name for Phot. En. / Polar. Rate waves"
prompt nmFolderMeas,"Name of the Folder with Mag. Meas. data"
prompt nmCoreNamesMeas,"Wave with Core Names of Mag. Meas. waves"
//prompt sufNamesMeas,"Name Suffix of Mag. Meas. waves"
prompt nmElec,"Electron Beam structure",popup wavelist("*"+SrwElecType,";","")
//prompt nmObs,"Radiation Sampling structure",popup wavelist("*"+SrwSmpType,";","")
prompt nmObsList,"List of Radiation Sampling structures " //,popup wavelist("*"+SrwSmpType,";","")
prompt numHarm,"Number of Harmonics to treat", popup "1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;10;21;22;23;24;25;26;27;28;29;30;31"
prompt relMinPhEn,"Rel. Start Fund. Phot. Energy Search"
prompt relMaxPhEn,"Rel. End Fund. Phot. Energy Search"
//prompt intMeth,SrwPMode,popup SrwPOPUPMode
prompt intPrec,"Relative Prec. for Integration along trajectory"
prompt rangePer,"Approx. Length of Periodic Field [m]"
Silent 1						|	Computing Radiation  ...
//PauseUpdate

srwUtiSetValS("nmCoreRes",nmCoreRes,"SrwUndPhotEnPolarVsParam")
srwUtiSetValS("nmFolderMeas",nmFolderMeas,"SrwUndPhotEnPolarVsParam")
srwUtiSetValS("nmCoreNamesMeas",nmCoreNamesMeas,"SrwUndPhotEnPolarVsParam")
srwUtiSetValS("nmObsList",nmObsList,"SrwUndPhotEnPolarVsParam")
//srwUtiSetValS("sufNamesMeas",sufNamesMeas,"SrwUndPhotEnPolarVsParam")
srwUtiSetValN("numHarm",numHarm,"SrwUndPhotEnPolarVsParam")
srwUtiSetValN("relMinPhEn",relMinPhEn,"SrwUndPhotEnPolarVsParam")
srwUtiSetValN("relMaxPhEn",relMaxPhEn,"SrwUndPhotEnPolarVsParam")
//srwUtiSetValN("intMeth",intMeth,"SrwUndPhotEnPolarVsParam")
srwUtiSetValN("intPrec",intPrec,"SrwUndPhotEnPolarVsParam")
srwUtiSetValN("rangePer",rangePer,"SrwUndPhotEnPolarVsParam")

if(strlen(nmCoreRes) > 26)
	abort "Core name for the resulting waves is too long"
endif

variable intMeth = 2 //"auto undulator"
variable precPer = 1.5*(0.003/intPrec)
if(precPer < 1.5)
	precPer = 1.5
endif

variable numObsStructs = dimsize($nmObsList, 0)
if(numObsStructs < 1)
	abort "Observation structures are not defined"
endif
string nmObsPinhole = $nmObsList[0]
string nmObsFinAp = ""
if(numObsStructs > 1)
	nmObsFinAp = $nmObsList[1]
endif
string nmObs = nmObsPinhole

string nmResPhotEnVsParam = nmCoreRes + "PhEn"
string nmResPolRateVsParam = nmCoreRes + "PolR"
string nmResPhotEnCorPin = nmCoreRes + "PhEnCpin"
string nmResPhotEnCorAp = nmCoreRes + "PhEnCap"

variable numParam1Vals = dimsize($nmCoreNamesMeas, 0)
variable numParam2Vals = dimsize($nmCoreNamesMeas, 1)
if(numParam2Vals <= 0)
	numParam2Vals = 1
endif

string strTest = nmFolderMeas[strlen(nmFolderMeas) - 1]
if(cmpstr(strTest, ":") != 0)
	nmFolderMeas += ":"
endif
string curNameBX, curNameBZ, curNameCore
variable sStart, sStep, ns, sRange, sOfstCp, halfRangeS

string AuxFldName = "AuxUndFundPhotEnPolar", AuxFldPerName = "AuxUndPhotEnPolPer"
string AuxFldNameBX = AuxFldName + "BX_fld", AuxFldNameBZ = AuxFldName + "BZ_fld"
string AuxFldNameTot = AuxFldName + SrwFieldType, AuxFldPerNameTot = AuxFldPerName + SrwUndType
string AuxShortFldName = "AuxUndFundPhotEnPolarS"
string AuxShortFldNameBX = AuxShortFldName + "BX_fld", AuxShortFldNameBZ = AuxShortFldName + "BZ_fld"
string AuxShortFldNameTot = AuxShortFldName + SrwFieldType
string AuxWfrName = "AuxUndFundPhotEnPol"
string AuxWfrNameTot = AuxWfrName + SrwRadType
string AuxSufInt = "Itot", AuxSufRlin = "Rlin", AuxSufRcirc = "Rcirc"
string AuxIntName = AuxWfrName + AuxSufInt + "_e"
string AuxRlinName = AuxWfrName + AuxSufRlin + "_e"
string AuxRcircName = AuxWfrName + AuxSufRcirc + "_e"
string AuxWaveToSearchFund = "AuxWaveToSearchFundPhtEn"
string AuxFldIntName = "AuxUndFundPhotEnFldInt"
string AuxDriftElecNameCore = "AuxDriftElec"
string AuxDriftElecNameTot = AuxDriftElecNameCore + "_bli"
string AuxStokesFilam = "AuxStokesFilam"
string AuxStokesFilamTot = AuxStokesFilam + "_ras"
string AuxStokesFilamIntSuf = "I0"
string AuxStokesFilamIntSpec = AuxStokesFilam + AuxStokesFilamIntSuf + "_e"
string AuxStokesFinEmPin = "AuxStokesFinEmPin"
string AuxStokesFinEmPinTot = AuxStokesFinEmPin + "_ras"
string AuxStokesFinEmPinIntSuf = "I1"
string AuxStokesFinEmPinIntSpec = AuxStokesFinEmPin + AuxStokesFinEmPinIntSuf + "_e"
string AuxStokesFinEmAp = "AuxStokesFinEmAp"
string AuxStokesFinEmApTot = AuxStokesFinEmAp + "_ras"
string AuxStokesFinEmApIntSuf = "I2"
string AuxStokesFinEmApIntSpec = AuxStokesFinEmAp + AuxStokesFinEmApIntSuf + "_e"

//variable factInfFundPhotEn = 2 // to tune: factor of the inferior limit, which determines the superior limit of the fundamental photon energy search
variable avgBXpeak, avgBZpeak, effUndPer, infFundPhotEn, supFundPhotEn, approxResonPhotEn, approxResonEn, twoFundPhotEn

variable obsPhotEnNp = srwGetSmpPhotEnNp(nmObs)
variable obsPhotEnStart = srwGetSmpPhotEnStart(nmObs)
variable obsPhotEnEnd = srwGetSmpPhotEnEnd(nmObs)
variable obsPhotEnStep = (obsPhotEnEnd - obsPhotEnStart)/(obsPhotEnNp - 1)
variable searchPhotEnNp

variable elecEn = srwGetElecBeamEnergy(nmElec)
variable elecS0 = srwGetElecBeamLongPos(nmElec)

variable auxI1X, auxI2X, auxI1Z, auxI2Z
variable kickSigmaS_m = 0.1 //2*approxUndPer*0.001 //to tune
variable dist_bw_Kicks, kickEntryHor, kickEntryVert, kickExitHor, kickExitVert
variable numHarmToCalc = (numHarm - 1)/2 + 1 //odd harmonics only

if(numParam2Vals == 1)
	make/O/N=(numParam1Vals, numHarmToCalc) $nmResPhotEnVsParam
	make/O/N=(numParam1Vals, 2, numHarmToCalc) $nmResPolRateVsParam //linear and circ.
else
	make/O/N=(numParam1Vals, numParam2Vals, numHarmToCalc) $nmResPhotEnVsParam
	make/O/N=(numParam1Vals, numParam2Vals, 2, numHarmToCalc) $nmResPolRateVsParam //linear and circ.
endif

duplicate/O $nmResPhotEnVsParam $nmResPhotEnCorPin
if(numObsStructs > 1)
	duplicate/O $nmResPhotEnVsParam $nmResPhotEnCorAp
endif

//make/O/N=(numParam1Vals, numParam2Vals, numHarmToCalc) $nmResPhotEnVsParam
//make/O/N=(numParam1Vals, numParam2Vals, 2, numHarmToCalc) $nmResPolRateVsParam //linear and circ.
//make/O/N=(numParam1Vals, numParam2Vals) $nmResPhotEnVsParam
//make/O/N=(numParam1Vals, numParam2Vals, 2) $nmResPolRateVsParam //linear and circ.

variable BxAndBzAreSupplied = 0
if(dimsize($nmCoreNamesMeas, 2) == 2)
	BxAndBzAreSupplied = 1
endif
if(BxAndBzAreSupplied == 0) //assuming only Bz is defined; creating zero Bx wave
	curNameBZ = nmFolderMeas + $nmCoreNamesMeas[0][0] //+ sufNamesMeas
	curNameBX = "AuxZeroUndPhotEnPolBx"
	duplicate/O $curNameBZ $curNameBX
	$curNameBX = 0
endif

string nmElecOrig = nmElec
nmElec = "Aux" + nmElecOrig
string nmElecFilam =  "AuxFil" + nmElecOrig
duplicate/O $nmElecOrig $nmElecFilam
SrwElecThick(nmElecFilam,1e-07,1e-07,1e-07,1,1,0,0,0,0)

variable npAuxFld = 20000
variable elecDriftLen
variable d_infFundPhotEn, d_supFundPhotEn
variable iParam2 = 0, iParam1, iHarm, indHarm
variable ePeakFilam, ePeakPin, ePeakAp
variable avgUndPer = 0, maxUndPer = 0, passCount = 0
do
	iParam1 = 0
	do
		//curNameCore = $nmCoreNamesMeas[iParam1][iParam2]+ sufNamesMeas
		//curNameBX = nmFolderMeas + "BX" + curNameCore
		//curNameBZ = nmFolderMeas + "BZ" + curNameCore
	
		if(BxAndBzAreSupplied == 0)
			curNameBZ = nmFolderMeas + $nmCoreNamesMeas[iParam1][iParam2] //+ sufNamesMeas
		else
			curNameBX = nmFolderMeas + $nmCoreNamesMeas[iParam1][iParam2][0] //+ sufNamesMeas
			curNameBZ = nmFolderMeas + $nmCoreNamesMeas[iParam1][iParam2][1] //+ sufNamesMeas
		endif
	
		sStart = dimoffset($curNameBZ, 0)
		sStep = dimdelta($curNameBZ, 0)
		ns = dimsize($curNameBZ, 0)
		sRange = sStep*(ns - 1)
		halfRangeS = 0.5*sRange
		sOfstCp = halfRangeS + sStart
		
		SrwMagFieldCreate(AuxFldName, 0, sRange, npAuxFld)
		$AuxFldNameBX = $curNameBX(x + sOfstCp)
		$AuxFldNameBZ = $curNameBZ(x + sOfstCp)
		
		//if(elecS0 < -halfRangeS)
		//	srwSetElecBeamLongPos(nmElec, -halfRangeS)
		//endif
		//if(elecS0 > halfRangeS)
		//	srwSetElecBeamLongPos(nmElec, halfRangeS)
		//endif
		
		elecDriftLen = -halfRangeS - elecS0
		SrwOptDrift(AuxDriftElecNameCore, elecDriftLen)
		duplicate/O $nmElecOrig $nmElec
		srElecBeamPropag($nmElec, $AuxDriftElecNameTot)
				
		//correcting field integrals to ensure right location of the axis of emission
		duplicate/O $AuxFldNameBX $AuxFldIntName
		integrate/T $AuxFldIntName
		auxI1X = $AuxFldIntName[npAuxFld - 1] //[T.m]
		integrate/T $AuxFldIntName
		auxI2X = $AuxFldIntName[npAuxFld - 1] //[T.m^2]
		
		duplicate/O $AuxFldNameBZ $AuxFldIntName
		integrate/T $AuxFldIntName
		auxI1Z = $AuxFldIntName[npAuxFld - 1] //[T.m]
		integrate/T $AuxFldIntName
		auxI2Z = $AuxFldIntName[npAuxFld - 1] //[T.m^2]

		dist_bw_Kicks = sRange - 6*kickSigmaS_m
		kickEntryHor = 0.5*(sRange/dist_bw_Kicks - 1)*auxI1Z - auxI2Z/dist_bw_Kicks
		kickExitHor = -0.5*(sRange/dist_bw_Kicks + 1)*auxI1Z + auxI2Z/dist_bw_Kicks
		kickEntryVert = 0.5*(sRange/dist_bw_Kicks - 1)*auxI1X - auxI2X/dist_bw_Kicks
		kickExitVert = -0.5*(sRange/dist_bw_Kicks + 1)*auxI1X + auxI2X/dist_bw_Kicks

		SrwMagGsnAng(AuxFldNameBX, 2,  -halfRangeS + 3*kickSigmaS_m, kickSigmaS_m*1000, kickEntryVert*1000)
		SrwMagGsnAng(AuxFldNameBX, 2,  halfRangeS - 3*kickSigmaS_m, kickSigmaS_m*1000, kickExitVert*1000)
		SrwMagGsnAng(AuxFldNameBZ, 2,  -halfRangeS + 3*kickSigmaS_m, kickSigmaS_m*1000, kickEntryHor*1000)
		SrwMagGsnAng(AuxFldNameBZ, 2,  halfRangeS - 3*kickSigmaS_m, kickSigmaS_m*1000, kickExitHor*1000)

		//calculating on-axis single-electron spectrum
		SrwMagPrec(AuxFldNameTot, intMeth, intPrec, intPrec, 10000, 1, 0, 0)
		SrwWfrCreate(AuxWfrName, nmElec, AuxFldNameTot, nmObs, 1, 1)
		
		SrwWfr2Int(AuxWfrNameTot, AuxSufInt, 7, 1, 1, 1, 1, 0, 0, 1) //on-axis total single-e intensity vs photon energy
		SrwWfr2PolRateExt(AuxWfrNameTot, AuxSufRlin, 1, 1, 1, 2, 1, 0, 0, 1) //on-axis single-e linear (hor./vert.) polarization rate vs photon energy
		SrwWfr2PolRateExt(AuxWfrNameTot, AuxSufRcirc, 5, 1, 1, 2, 1, 0, 0, 1) //on-axis single-e circular (right/left) polarization rate vs photon energy
		
		maxUndPer = 0.25*sRange*1000
		if(avgUndPer != 0)
			maxUndPer = 1.2*avgUndPer
		endif
		
		if(rangePer > 0)
			SrwMagDupl(AuxFldNameTot, AuxShortFldName)
			$AuxShortFldNameBX = 0; $AuxShortFldNameBZ = 0
			$AuxShortFldNameBX = $AuxFldNameBX(x)*srwUtiNonZeroInterval(x, -0.5*rangePer, 0.5*rangePer)
			$AuxShortFldNameBZ = $AuxFldNameBZ(x)*srwUtiNonZeroInterval(x, -0.5*rangePer, 0.5*rangePer)
			SrwMagArb2Per(AuxFldPerName, AuxShortFldNameTot, 0.05, 3, maxUndPer)
		else
			SrwMagArb2Per(AuxFldPerName, AuxFldNameTot, 0.05, 3, maxUndPer)
		endif
				
		effUndPer = str2num($AuxFldPerNameTot[0]) //undulator period in [m]
		avgUndPer = (passCount*avgUndPer + effUndPer)/(passCount + 1)

		//estimating value of the fundamental photon energy
		//to improve: get Kx, Kz from periodic field structure
		wavestats/Q/R=(-2*avgUndPer, 2*avgUndPer) $AuxFldNameBX
		avgBXpeak = 0.5*(abs(V_min) + abs(V_max))
		wavestats/Q/R=(-2*avgUndPer, 2*avgUndPer) $AuxFldNameBZ
		avgBZpeak = 0.5*(abs(V_min) + abs(V_max))

		SrwPerStoCreate(AuxStokesFilam,nmElecFilam,AuxFldPerNameTot,nmObsPinhole,1,numHarm+2,precPer,precPer,2)
		SrwSto2IntF(AuxStokesFilamTot,AuxStokesFilamIntSuf,7,1,1,1,0,0,1)

		SrwPerStoCreate(AuxStokesFinEmPin,nmElecOrig,AuxFldPerNameTot,nmObsPinhole,1,numHarm+2,precPer,precPer,2)
		SrwSto2IntF(AuxStokesFinEmPinTot,AuxStokesFinEmPinIntSuf,7,1,1,1,0,0,1)
		
		if(strlen(nmObsFinAp) > 0)
			SrwPerStoCreate(AuxStokesFinEmAp,nmElecOrig,AuxFldPerNameTot,nmObsFinAp,1,numHarm+2,precPer,precPer,1)
			SrwSto2IntF(AuxStokesFinEmApTot,AuxStokesFinEmApIntSuf,7,1,1,1,0,0,1)
		endif
		
		//infFundPhotEn = srUtiUndFundPhotEn(sqrt(avgBXpeak*avgBXpeak + avgBZpeak*avgBZpeak), effUndPer, elecEn, 2) // [eV]
		//supFundPhotEn = infFundPhotEn*factInfFundPhotEn
		
		approxResonPhotEn = srUtiUndFundPhotEn(sqrt(avgBXpeak*avgBXpeak + avgBZpeak*avgBZpeak), effUndPer, elecEn, 2) // [eV]
		infFundPhotEn = relMinPhEn*approxResonPhotEn
		supFundPhotEn = relMaxPhEn*approxResonPhotEn
		searchPhotEnNp = round((supFundPhotEn - infFundPhotEn)/obsPhotEnStep + 1)
		make/O/N=(searchPhotEnNp) $AuxWaveToSearchFund
		
		d_infFundPhotEn  = approxResonPhotEn - infFundPhotEn
		d_supFundPhotEn = supFundPhotEn - approxResonPhotEn
		approxResonEn = approxResonPhotEn
		
		iHarm = 1
		do
			infFundPhotEn = approxResonEn - d_infFundPhotEn
			supFundPhotEn = approxResonEn + d_supFundPhotEn
			
			if(infFundPhotEn > obsPhotEnEnd)
				print "WARNING: Estimated resonant photon energy (" , approxResonEn, "eV)  at harmonic", iHarm, "is out of range"
				break
			endif
			if(supFundPhotEn < obsPhotEnStart)
				print "WARNING: Estimated resonant photon energy (" , approxResonEn, "eV)  at harmonic", iHarm, "is out of range"
			endif
			
			SetScale/I x infFundPhotEn, supFundPhotEn, "eV", $AuxWaveToSearchFund
			$AuxWaveToSearchFund = $AuxIntName(x)
			wavestats/Q $AuxWaveToSearchFund
			
			if(V_maxloc <= infFundPhotEn)
				print "WARNING: Lower limit reached when searching harmonic No:", iHarm
			endif
			if(V_maxloc >= supFundPhotEn)
				print "WARNING: Upper limit reached when searching harmonic No:", iHarm
			endif
			
			if(V_maxloc >= obsPhotEnEnd)
				print "WARNING: Upper limit of observation photon energy was reached when searching harmonic No:", iHarm
				break
			endif
			
			indHarm = round((iHarm - 1)/2)
			
			if(numParam2Vals == 1)
				$nmResPhotEnVsParam[iParam1][indHarm] = V_maxloc //fundamental photon energy
				$nmResPolRateVsParam[iParam1][0][indHarm] = $AuxRlinName(V_maxloc) //linear polarization rate on fundamental
				$nmResPolRateVsParam[iParam1][1][indHarm] = $AuxRcircName(V_maxloc) //circular polarization rate on fundamental
			else
				$nmResPhotEnVsParam[iParam1][iParam2][indHarm] = V_maxloc //fundamental photon energy
				$nmResPolRateVsParam[iParam1][iParam2][0][indHarm] = $AuxRlinName(V_maxloc) //linear polarization rate on fundamental
				$nmResPolRateVsParam[iParam1][iParam2][1][indHarm] = $AuxRcircName(V_maxloc) //circular polarization rate on fundamental
			endif
			
			if(iHarm == 1)
				twoFundPhotEn = V_maxloc*2
			endif
			approxResonEn = V_maxloc + twoFundPhotEn
			
			$AuxWaveToSearchFund = $AuxStokesFilamIntSpec(x)
			wavestats/Q $AuxWaveToSearchFund
			ePeakFilam = V_maxloc
			
			$AuxWaveToSearchFund = $AuxStokesFinEmPinIntSpec(x)
			wavestats/Q $AuxWaveToSearchFund
			ePeakPin = V_maxloc
			
			if(numParam2Vals == 1)
				$nmResPhotEnCorPin[iParam1][indHarm] = ePeakFilam - ePeakPin
			else
				$nmResPhotEnCorPin[iParam1][iParam2][indHarm] = ePeakFilam - ePeakPin
			endif
			
			if(strlen(nmObsFinAp) > 0)
				$AuxWaveToSearchFund = $AuxStokesFinEmApIntSpec(x)
				wavestats/Q $AuxWaveToSearchFund
				ePeakAp = V_maxloc
				
				if(numParam2Vals == 1)
					$nmResPhotEnCorAp[iParam1][indHarm] = ePeakFilam - ePeakAp
				else
					$nmResPhotEnCorAp[iParam1][iParam2][indHarm] = ePeakFilam - ePeakAp
				endif
			endif
						
			iHarm += 2
		while(iHarm <= numHarm)

			//DoUpdate
		passCount += 1
		iParam1 += 1
	while(iParam1 < numParam1Vals)
	iParam2 += 1
while(iParam2 < numParam2Vals)

srwUtiSetValS("SrwElecName",nmElec[0,strlen(nmElec)-strlen(SrwElecType)-1],"")
srwUtiSetValS("SrwSmpName",nmObs[0,strlen(nmObs)-strlen(SrwSmpType)-1],"")
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Generating final operation tables: Gap vs Photon Energy
//Usage: SrwUndGapVsPhotEn("otGapVsPhEn_U20_CRISTAL", "otU20_CRISTALPhEn", "GapVals_U20_CRISTAL", 1400, 100, 237)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUndGapVsPhotEn(nmResCore, nmDepPhotEnVsUndGap, nmUndGapVals, eStart, eStep, eNp)
string nmResCore, nmDepPhotEnVsUndGap, nmUndGapVals
variable eStart, eStep, eNp

variable numUndGapVals = dimsize($nmDepPhotEnVsUndGap, 0)
variable numUndParam2Vals = dimsize($nmDepPhotEnVsUndGap, 1)
variable numHarmInd = dimsize($nmDepPhotEnVsUndGap, 2)
if(numHarmInd <= 0)
	numHarmInd = dimsize($nmDepPhotEnVsUndGap, 1)
	numUndParam2Vals = 1
endif

make/D/N=5/O W_coef
string nmResHarm
variable indHarm = 0, gapMin = $nmUndGapVals[0], gapMax = $nmUndGapVals[dimsize($nmUndGapVals, 0) - 1]
variable numActGapVals
do
	make/O/N=(numUndGapVals) wAuxPhotEnVsUndGapAtHarm
	duplicate/O $nmUndGapVals wAuxGapValsOrig
	
	if(numUndParam2Vals <= 1)
		wAuxPhotEnVsUndGapAtHarm = $nmDepPhotEnVsUndGap[p][indHarm]
	else
	
	endif
	
	numActGapVals = 1
	do
		if(wAuxPhotEnVsUndGapAtHarm[numActGapVals - 1] <= 0)
			numActGapVals -= 1
			break
		endif
		numActGapVals += 1
	while(numActGapVals <= numUndGapVals)
	
	if(numActGapVals < numUndGapVals)
		redimension/N=(numActGapVals) wAuxPhotEnVsUndGapAtHarm
		redimension/N=(numActGapVals) wAuxGapValsOrig
	endif

	nmResHarm = nmResCore + num2str(2*indHarm + 1)
	make/O/N=(eNp) $nmResHarm
	SetScale/P x eStart,eStep,"", $nmResHarm

	if(numActGapVals >= 5)
		W_coef[0] = {0,wAuxPhotEnVsUndGapAtHarm[10000],1,0,1}
		FuncFit srwUndPhotEnVsGap W_coef  wAuxPhotEnVsUndGapAtHarm /X=wAuxGapValsOrig /D 
		$nmResHarm = srwUndFindGapForPhotEn(W_coef, x, gapMin, gapMax)
	else
		if(numActGapVals >= 2)
			Interpolate2/T=1/N=1000/Y=wAuxPhotEnVsUndGapAtHarm_L wAuxGapValsOrig, wAuxPhotEnVsUndGapAtHarm
			$nmResHarm = srwUndFindGapForPhotEnInterp(wAuxPhotEnVsUndGapAtHarm_L, x, gapMin, gapMax)
		endif
	endif

	indHarm += 1
while(indHarm < numHarmInd)

end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Function fitting dependence of Phot. Energy vs Gap (Sigmoid)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function srwUndPhotEnVsGap(w,x) : FitFunc
	Wave w
	Variable x

	//return w[0] + w[1]/(1 +(1+ w[4]*x)*exp(-(x - w[2])/w[3]))
	//return w[0] + w[1]/(1 + exp(-((x - w[2])*(1 + w[4]*(x - w[2])))/w[3]))
	//return w[0] + w[1]/(1 + exp(-(x - w[2])/w[3]))
	
	return w[0] + w[1]/(1 + w[2]*exp(-(x - w[3])/w[4]))
	
	//variable/G g0_PhotEnVsGap, e0_PhotEnVsGap
	//variable w0 = e0_PhotEnVsGap - w[0]/(1 + exp(-(x - w[1])/w[2]))
End

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Function fitting dependence of Phot. Energy vs Phase (Polygon)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function srwUndPhotEnVsPhase(w,x) : FitFunc
	wave w
	variable x
	return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*(w[4]+x*(w[5]+x*(w[6]+x*(w[7]+x*w[8])))))))

	//return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*(w[4]+x*(w[5]+x*(w[6]+x*(w[7]+x*(w[8]+x*(w[9]+x*w[10])))))))))
	//return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*(w[4]+x*w[5]))))
	//return w[0] + x*(w[1] +  x*(w[2] +x*w[3]))

End

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Function fitting dependence of Polarization Rate vs Gap (Poly)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function srwUndPolRateVsGap(w,x) : FitFunc
	wave w
	variable x
	//return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*w[4])))
	return w[0] + x*(w[1] +  x*(w[2] +x*w[3]))
End

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Function fitting dependence of Polarization Rate vs Phase (Poly)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function srwUndPolRateVsPhase(w,x) : FitFunc
	wave w
	variable x
	//return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*(w[4] + x*(w[5] + x*(w[6] + x*(w[7]+x*w[8])))))))
	return w[0] + x*(w[1] +  x*(w[2] +x*(w[3] + x*(w[4] + x*(w[5] + x*w[6])))))
End

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Finds gap for given photon energy
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndFindGapForPhotEn(w, PhotEn, GapMin, GapMax)
wave w
variable PhotEn, GapMin, GapMax
//variable/G V_Root = 0 //OC commented out at porting to Igor6
FindRoots/Q /L=(GapMin) /H=(GapMax) /Z=(PhotEn) srwUndPhotEnVsGap, w

if(V_Root < GapMin)
	V_Root = 0
endif
if(V_Root > GapMax)
	V_Root = 0
endif
return V_Root
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Finds gap for given photon energy by linear interpolation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndFindGapForPhotEnInterp(w, PhotEn, GapMin, GapMax)
wave w
variable PhotEn, GapMin, GapMax
variable gapStart = dimoffset(w, 0)
variable gapEnd = gapStart + (dimsize(w, 0) - 1)*dimdelta(w, 0)
//variable/G V_LevelX = NaN //OC commented out at porting to Igor6
FindLevel /Q/R=(GapMin, GapMax) w, PhotEn

if((V_LevelX >= gapStart) %& (V_LevelX <= gapEnd))
	return V_LevelX
else
	return 0
endif

//if(V_LevelX == NaN)
//	return 0
//else
//	if((V_LevelX < gapStart) %| (V_LevelX > gapEnd))
//		return 0
//	endif
//endif
//return V_LevelX
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndPhotEnVsGapAndPhase(wCoefsGapPhase, gap, phase)
wave wCoefsGapPhase
variable gap, phase

variable nCoefsGap = dimsize(wCoefsGapPhase, 0)
variable nCoefsPhase = dimsize(wCoefsGapPhase, 1)
make/O/N=(nCoefsGap) wAuxCoefsGap
make/O/N=(nCoefsPhase) wAuxCoefsPhase

variable i=0
do
	wAuxCoefsPhase = wCoefsGapPhase[i][p]
	wAuxCoefsGap[i] = srwUndPhotEnVsPhase(wAuxCoefsPhase, phase)
	i += 1
while(i < nCoefsGap)
return srwUndPhotEnVsGap(wAuxCoefsGap, gap)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndPolRateVsGapAndPhase(wCoefsGapPhase, gap, phase)
wave wCoefsGapPhase
variable gap, phase

variable nCoefsGap = dimsize(wCoefsGapPhase, 0)
variable nCoefsPhase = dimsize(wCoefsGapPhase, 1)
make/O/N=(nCoefsGap) wAuxCoefsGap
make/O/N=(nCoefsPhase) wAuxCoefsPhase

variable i=0
do
	wAuxCoefsPhase = wCoefsGapPhase[i][p]
	wAuxCoefsGap[i] = srwUndPolRateVsPhase(wAuxCoefsPhase, phase)
	i += 1
while(i < nCoefsGap)
return srwUndPolRateVsGap(wAuxCoefsGap, gap)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndDevPhotEnAndPolRate(wPhotEnAndPolRate, gap, phase)
wave wPhotEnAndPolRate
variable gap, phase

variable goalPhotEn = wPhotEnAndPolRate[0]
variable goalPolRate = wPhotEnAndPolRate[1]

//These waves should EXIST !!! Create them before calling this function !!!
wave otPhotEnTestCoefsAux //to make more general
wave otPolRateTestCoefsAux

variable dPhotEn = srwUndPhotEnVsGapAndPhase(otPhotEnTestCoefsAux, gap, phase) - goalPhotEn
variable dPolRate = srwUndPolRateVsGapAndPhase(otPolRateTestCoefsAux, gap, phase) - goalPolRate

//return sqrt(dPhotEn*dPhotEn/(goalPhotEn*goalPhotEn) + dPolRate*dPolRate/(goalPolRate*goalPolRate))
//return sqrt(dPhotEn*dPhotEn + 5*dPolRate*dPolRate)
//return sqrt(dPhotEn*dPhotEn + 5*dPolRate*dPolRate)
return dPhotEn*dPhotEn + 10*dPolRate*dPolRate

end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUndDevPhotEnVsGapForPhase(wPhotEn, gap)
wave wPhotEn
variable gap

variable goalPhotEn = wPhotEn[0]
//variable goalPolRate = wPhotEnAndPolRate[1]

//These waves should EXIST !!! Create them before calling this function !!!
wave otPhotEnTestCoefsAux //to make more general
wave otPhaseCoefsAux
variable polyOrder = dimsize(otPhaseCoefsAux, 0)
variable phaseVal =  poly(otPhaseCoefsAux, goalPhotEn)

variable dPhotEn = srwUndPhotEnVsGapAndPhase(otPhotEnTestCoefsAux, gap, phaseVal) - goalPhotEn
//variable dPolRate = srwUndPolRateVsGapAndPhase(otPolRateTestCoefsAux, gap, phase) - goalPolRate

//return sqrt(dPhotEn*dPhotEn/(goalPhotEn*goalPhotEn) + dPolRate*dPolRate/(goalPolRate*goalPolRate))
//return sqrt(dPhotEn*dPhotEn + 5*dPolRate*dPolRate)
//return sqrt(dPhotEn*dPhotEn + 5*dPolRate*dPolRate)
return dPhotEn*dPhotEn
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUndFindGapAndPhaseFromParam(nmResGap, nmResPhase, nmCoefPhotEn, nmCoefPolRate, nmPhotEn, polRate, nmInitGapPhase)
string nmResGap, nmResPhase
string nmCoefPhotEn, nmCoefPolRate
string nmPhotEn, nmInitGapPhase
variable polRate
Silent 1						|	Computing Radiation  ...

duplicate/O $nmPhotEn $nmResGap; $nmResGap = 0
duplicate/O $nmPhotEn $nmResPhase; $nmResPhase = 0

variable numPhotEnVals = dimsize($nmPhotEn, 0)
//make/O/D wGapPhaseFound = {15., 15}
duplicate/O $nmInitGapPhase wGapPhaseFound
make/O/N=2 wPhotEnAndPolar
wPhotEnAndPolar = polRate

duplicate/O $nmCoefPhotEn otPhotEnTestCoefsAux //because the function srwUndDevPhotEnAndPolRate needs these waves
duplicate/O $nmCoefPolRate otPolRateTestCoefsAux

variable i=0
do
	wPhotEnAndPolar[0] = $nmPhotEn[i]
	Optimize/Q /T={8.53618E-09, 7.28664E-14}/D=15 /X=wGapPhaseFound srwUndDevPhotEnAndPolRate, wPhotEnAndPolar
	//Optimize/Q /T={8.53618E-09, 7.28664E-14} /M = {1, 1 } /X=wGapPhaseFound srwUndDevPhotEnAndPolRate, wPhotEnAndPolar
	//Optimize/Q /X=wGapPhaseFound srwUndDevPhotEnAndPolRate, wPhotEnAndPolar
	$nmResGap[i] = wGapPhaseFound[0]
	$nmResPhase[i] = wGapPhaseFound[1]
	
		DoUpdate

	i += 1
while(i < numPhotEnVals)

//killwaves/Z otPhotEnTestCoefsAux, otPolRateTestCoefsAux
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary; to make more general
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUndFindGapForGivenPhase(nmResGap, nmCoefPhotEn, nmCoefPhase, nmPhotEn, gapMin, gapMax)
string nmResGap //, nmResPhase
string nmCoefPhotEn, nmCoefPhase
string nmPhotEn
variable gapMin, gapMax
Silent 1						|	Computing Radiation  ...

duplicate/O $nmPhotEn $nmResGap; $nmResGap = 0

variable numPhotEnVals = dimsize($nmPhotEn, 0)
//make/O/D wGapFound = {15.5}, wPhotEn
make/O/D wPhotEn

duplicate/O $nmCoefPhotEn otPhotEnTestCoefsAux //because the function srwUndDevPhotEnAndPolRate needs these waves
duplicate/O $nmCoefPhase otPhaseCoefsAux

variable i=0
do
	wPhotEn[0] = $nmPhotEn[i]
	//Optimize/Q /T={8.53618E-09, 7.28664E-14}/D=15 /X=wGapFound srwUndDevPhotEnVsGapForPhase, wPhotEn
	Optimize/Q /T=0.00001/D=15 /L=(gapMin) /H=(gapMax) srwUndDevPhotEnVsGapForPhase, wPhotEn

	//$nmResGap[i] = wGapFound[0]
	$nmResGap[i] = V_minloc
	
		DoUpdate

	i += 1
while(i < numPhotEnVals)

killwaves/Z otPhotEnTestCoefsAux, otPhaseCoefsAux
end

