
//+++++++++++++++++++++++++++++++++++++++
//
//Propagate wavefront to deduce Stokes params of radiation emitted
//by Thick electron beam.
//In this version, initial wavefront is re-calculated for every macro-particle.
//Propagation through each optical element of a Beamline can be controlled.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrEmitPropStokesMultiE(StoName, ElecName, MagName, ObsName, ObsNxNzSamplFact, PropBL, MaxPrt, radCmpn2View, xcView, zcView)
string StoName=srwUtiGetValS("StoName", "Stk", "SrwWfrEmitPropStokesMultiE")
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagName+SrwFieldType
string ObsName=SrwSmpName+SrwSmpType
variable ObsNxNzSamplFact=srwUtiGetValN("ObsNxNzSamplFact", SrwSmpNxNzSamplFact, "SrwWfrEmitPropStokesMultiE") 
string PropBL=srwUtiGetValS("PropBL", "", "SrwWfrEmitPropStokesMultiE") //SrwBliLast+SrwBeamlineType
variable MaxPrt=srwUtiGetValN("MaxPrt", 10000, "SrwWfrEmitPropStokesMultiE")
variable radCmpn2View=srwUtiGetValN("radCmpn2View", 1, "SrwWfrEmitPropStokesMultiE")
variable xcView=srwUtiGetValN("xcView", 0, "SrwWfrEmitPropStokesMultiE")
variable zcView=srwUtiGetValN("zcView", 0, "SrwWfrEmitPropStokesMultiE")
prompt StoName, "Name of the Stokes structure"
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt ObsNxNzSamplFact,"Oversampling Factor (effective if > 0)"
prompt PropBL, "Optical Elem. and Propag. Param. List",popup WaveList("*",";","TEXT:1,DIMS:2,MINCOLS:10")
prompt MaxPrt, "Number of Macro-Particles"
prompt radCmpn2View, "Polarization Component to View", popup "-none-;"+SrwPOPUPPolar+";Total"
prompt xcView, "Hor. Center Point for Viewing [mm]"
prompt zcView, "Vert. Center Point for Viewing [mm]"
Silent 1						|	Propagating Wavefront ...
PauseUpdate

if(exists(PropBL) != 1)
	abort "List of Optical Elements and Propagation Parameters was not provided"
endif

SrwElecName = ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName = ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]

srwUtiSetValS("StoName", StoName, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValN("ObsNxNzSamplFact", ObsNxNzSamplFact, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValS("PropBL", PropBL, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValN("MaxPrt", MaxPrt, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValN("radCmpn2View", radCmpn2View, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValN("xcView", xcView, "SrwWfrEmitPropStokesMultiE")
srwUtiSetValN("zcView", zcView, "SrwWfrEmitPropStokesMultiE")

//-----important steering parameters!!!
string nmAccelDistrFirstAp = "wAccelDistr1stAp_xzE"
string nmAccelDistrWaist = "wAccelDistrWaist_xz"
string nmAccelDistrElecEn = "wFluxVsElecEnNorm"

//variable initWfrCalcMeth = 2 //3- "auto wiggler", 2- "auto undulator"
//variable initWfrCalcRelPrec = 0.01 //0.01
//variable relFluxRejectPropag = 1e-04 //rel. flux under which propagation should not be performed
//variable treatElecEnergySpreadOnly = 0 //0 by default //1

//variable angIntens = 0 //0- extract intensity in Coordinate representation, 1- in Angular representation
//variable azimAngForIntensCuts = 0 //Pi/4.
//variable logScaleIntens = 0 //1
//-----

//Eventually call another proc for more input parameters
variable justCalled = srwUtiGetValN("justCalled", 0, "SrwWfrEmitPropStokesAuxIn")
if(justCalled == 0)
	SrwWfrEmitPropStokesAuxIn()
endif
srwUtiSetValN("justCalled", 0, "SrwWfrEmitPropStokesAuxIn")
variable relFluxRejectPropag = srwUtiGetValN("relFluxRejectProp", 1e-04, "SrwWfrEmitPropStokesAuxIn") //rel. flux under which propagation should not be performed
variable treatElecEnergySpreadOnly = srwUtiGetValN("treatElecEnSpreadOnly", 1, "SrwWfrEmitPropStokesAuxIn") - 1  //0 by default //1
variable photEnIntegIntensUnit = srwUtiGetValN("photEnIntegIntensUnit", 1, "SrwWfrEmitPropStokesAuxIn")   //1- Ph/s/mm^2; 2- W/mm^2

variable angIntens = srwUtiGetValN("angIntens", 1, "SrwWfrEmitPropStokesAuxIn") - 1 //0- extract intensity in Coordinate representation, 1- in Angular representation
variable azimAngForIntensCuts = srwUtiGetValN("azimAngForIntensCuts", 0, "SrwWfrEmitPropStokesAuxIn")
variable logScaleIntens = srwUtiGetValN("scaleIntens", 1, "SrwWfrEmitPropStokesAuxIn") - 1 //0- use lusual inear scale, 1- use logarithmic scale

string nmSpecTr = srwUtiGetValS("nmSpecTr", "", "SrwWfrEmitPropStokesAuxIn")
if(exists(nmSpecTr) == 0)
	nmSpecTr = ""
endif
variable specTrIsDef = 1
if(strlen(nmSpecTr) <= 0)
	specTrIsDef = 0
endif

variable cosAzimAng = cos(azimAngForIntensCuts)
variable sinAzimAng = sin(azimAngForIntensCuts)

variable convPhW = 1.60218e-16 //1  Phot/s/.1%bw    correspond(s) to :   1.60218e-16  W/eV

variable AccelDistrFirstApExists = exists(nmAccelDistrFirstAp)
variable AccelDistrWaistExists = exists(nmAccelDistrWaist)
variable AccelDistrElecEnExists = exists(nmAccelDistrElecEn)

variable maxAccelDistrFirstAp = 1, maxAccelDistrWaist = 1, maxAccelDistrElecEn = 1
variable minElecEnAccelDistrFirstAp, maxElecEnAccelDistrFirstAp
variable minProjXAccelDistrFirstAp, maxProjXAccelDistrFirstAp
variable minProjZAccelDistrFirstAp, maxProjZAccelDistrFirstAp
variable minXAccelDistrWaist, maxXAccelDistrWaist
variable minZAccelDistrWaist, maxZAccelDistrWaist

if(AccelDistrFirstApExists)
	wavestats/Q $nmAccelDistrFirstAp
	maxAccelDistrFirstAp = V_max
	
	minProjXAccelDistrFirstAp = dimoffset($nmAccelDistrFirstAp, 0)
	maxProjXAccelDistrFirstAp = minProjXAccelDistrFirstAp + (dimsize($nmAccelDistrFirstAp, 0) - 1)*dimdelta($nmAccelDistrFirstAp, 0)
	minProjZAccelDistrFirstAp = dimoffset($nmAccelDistrFirstAp, 1)
	maxProjZAccelDistrFirstAp = minProjZAccelDistrFirstAp + (dimsize($nmAccelDistrFirstAp, 1) - 1)*dimdelta($nmAccelDistrFirstAp, 1)
	minElecEnAccelDistrFirstAp = dimoffset($nmAccelDistrFirstAp, 2)
	maxElecEnAccelDistrFirstAp = minElecEnAccelDistrFirstAp + (dimsize($nmAccelDistrFirstAp, 2) - 1)*dimdelta($nmAccelDistrFirstAp, 2)
endif
if(AccelDistrWaistExists)
	wavestats/Q $nmAccelDistrWaist
	maxAccelDistrWaist = V_max
	
	minXAccelDistrWaist = dimoffset($nmAccelDistrWaist, 0)
	maxXAccelDistrWaist = minXAccelDistrWaist + (dimsize($nmAccelDistrWaist, 0) - 1)*dimdelta($nmAccelDistrWaist, 0)
	minZAccelDistrWaist = dimoffset($nmAccelDistrWaist, 1)
	maxZAccelDistrWaist = minZAccelDistrWaist + (dimsize($nmAccelDistrWaist, 1) - 1)*dimdelta($nmAccelDistrWaist, 1)
endif
if(AccelDistrElecEnExists)
	wavestats/Q $nmAccelDistrElecEn
	maxAccelDistrElecEn = V_max
endif
variable halfMaxAccelDistrFirstAp = 0.5*maxAccelDistrFirstAp
variable halfMaxAccelDistrWaist = 0.5*maxAccelDistrWaist
variable halfMaxAccelDistrElecEn = 0.5*maxAccelDistrElecEn

variable ObsNxNzForProp = 1
if(ObsNxNzSamplFact > 0)
	ObsNxNzForProp = 2
endif

string origObsName = ObsName
string auxObsName = "AuxObsEmitPropStokesMultiE_obs"

variable startPhotEn = srwGetSmpPhotEnStart(ObsName) //[eV]
variable endPhotEn = srwGetSmpPhotEnEnd(ObsName) //[eV]
variable absPhotEnInterv = 0, halfAbsPhotEnInterv, cenPhotEn
if(startPhotEn != endPhotEn)
	absPhotEnInterv = endPhotEn - startPhotEn
	halfAbsPhotEnInterv = 0.5*absPhotEnInterv
	cenPhotEn = 0.5*(startPhotEn + endPhotEn)
	//intensity will integrated over photon energy within this absolute finite interval
	duplicate/O $ObsName $auxObsName
	ObsName = auxObsName
endif

variable elecSigX = srwGetElecBeamHorSizeRMS(ElecName)
variable elecSigXp = srwGetElecBeamHorDivergRMS(ElecName)
variable elecMXXp = srwGetElecBeamHorMixedMom(ElecName)
variable elecSigXe2 = elecSigX*elecSigX
variable elecSigXpe2 = elecSigXp*elecSigXp
variable elecSigZ = srwGetElecBeamVertSizeRMS(ElecName)
variable elecSigZp = srwGetElecBeamVertDivergRMS(ElecName)
variable elecMZZp = srwGetElecBeamVertMixedMom(ElecName)
variable elecSigZe2 = elecSigZ*elecSigZ
variable elecSigZpe2 = elecSigZp*elecSigZp
variable elecRelEnSpr = srwGetElecBeamRelEnSprRMS(ElecName)

variable elecEn0 = srwGetElecBeamEnergy(ElecName)
variable elecAbsEnSpr = elecEn0*elecRelEnSpr

variable elecCur = srwGetElecBeamCurrent(ElecName)
variable elecS0 = srwGetElecBeamLongPos(ElecName)
variable elecX0 = srwGetElecBeamHorPos(ElecName)
variable elecXp0 = srwGetElecBeamHorAng(ElecName)
variable elecZ0 = srwGetElecBeamVertPos(ElecName)
variable elecZp0 = srwGetElecBeamVertAng(ElecName)

variable elecLimitEstart = 0, elecLimitEend = 1000
variable elecLimitXstart = -1000, elecLimitXend = 1000
variable elecLimitXPstart = -1000, elecLimitXPend = 1000
variable elecLimitZstart = -1000, elecLimitZend = 1000
variable elecLimitZPstart = -1000, elecLimitZPend = 1000
variable hrE, hrX, hrXp, hrZ, hrZp
if(dimsize($ElecName, 0) >= 51)
	hrE = 0.5*$ElecName[46]
	if(hrE > 0)
		elecLimitEstart = elecEn0 - hrE
		elecLimitEend = elecEn0 + hrE
	endif
	hrX = 0.5*$ElecName[47]
	if(hrX > 0)
		elecLimitXstart = elecX0 - hrX
		elecLimitXend = elecX0 + hrX
	endif
	hrXp = 0.5*$ElecName[48]
	if(hrXp > 0)
		elecLimitXPstart = elecXp0 - hrXp
		elecLimitXPend = elecXp0 + hrXp
	endif
	hrZ = 0.5*$ElecName[49]
	if(hrZ > 0)
		elecLimitZstart = elecZ0 - hrZ
		elecLimitZend = elecZ0 + hrZ
	endif
	hrZp = 0.5*$ElecName[50]
	if(hrZp > 0)
		elecLimitZPstart = elecZp0 - hrZp
		elecLimitZPend = elecZp0 + hrZp
	endif
endif

string ElecWorkCore = "elecEmitPropStokesMultiE"
string ElecWorkName = ElecWorkCore + "_ebm"
string RadWorkCore = StoName //"wfrEmitPropStokesMultiE"
string RadWorkName = RadWorkCore + "_rad"
string nmRadWorkEX = RadWorkCore + "X_rae"
string nmRadWorkEZ = RadWorkCore + "Z_rae"
string StoWorkName = RadWorkCore + SrwStoType

string IntCurSuf = "Icur", IntResSuf = "Ires"
string nmXcSigLPTau = "wCentersSigmasPropStokesMultiE"
string nmElecInitCond = "wElecInitCondPropStokesMultiE"

string nmCurWorkIntXZ = RadWorkCore + IntCurSuf + "_xz"
string nmCurWorkIntX = RadWorkCore + IntCurSuf + "_x"
string nmCurWorkIntZ = RadWorkCore + IntCurSuf + "_z"

string nmCurStoWorkIntXZ = RadWorkCore + IntResSuf + "_xz"
string nmCurStoWorkIntX = RadWorkCore + IntResSuf + "_x"
string nmCurStoWorkIntZ = RadWorkCore + IntResSuf + "_z"

string sufRadWorkFlux = "F"
string nmRadWorkFlux = RadWorkCore + sufRadWorkFlux + "_e"

//SrwMagPrec(MagName,initWfrCalcMeth,0.01,initWfrCalcRelPrec,10000,1,0,0)
SrwMagPrec(MagName,SrwMode,SrwRadIntStep,SrwPrec,10000,1,0,0)
SrwElecFilament(ElecWorkCore, elecEn0, elecCur, elecS0, elecX0*1000, elecXp0*1000, elecZ0*1000, elecZp0*1000)

variable curPhotEn = startPhotEn
if(startPhotEn != endPhotEn)
	curPhotEn = cenPhotEn
	srwSetSmpPhotEnStart(ObsName, curPhotEn)
	srwSetSmpPhotEnEnd(ObsName, curPhotEn)
endif
SrwWfrCreate(RadWorkCore, ElecWorkName, MagName, ObsName, ObsNxNzForProp, abs(ObsNxNzSamplFact))

variable propagWasDone = 0

SrwWfr2Int(RadWorkName,sufRadWorkFlux,7,5,1,1,1,0,0,1) //flux
variable curFlux = 0, curMaxFlux = $nmRadWorkFlux[0], curPropFlux

SrwWfrPropList(RadWorkName, PropBL, 1, "")

if(angIntens)
	srWfrSetRepres($RadWorkName, "A")
endif

SrwWfr2Sto(RadWorkCore, RadWorkName)

variable curValAccelDistrFirstAp, testValAccelDistrFirstAp, elecProjX, elecProjZ, AccelDistrFirstApPasses
variable rProj = srwGetSmpLongPos(ObsName) - elecS0
variable curValAccelDistrWaist, testValAccelDistrWaist, AccelDistrWaistPasses, AccelDistrElecEnPasses
variable testValAccelDistrElecEn, curValAccelDistrElecEn
variable multDistr = 1

if(AccelDistrFirstApExists)
	elecProjX = -elecX0 - rProj*elecXp0
	elecProjZ = -elecZ0 - rProj*elecZp0
	curValAccelDistrFirstAp = $nmAccelDistrFirstAp(elecProjX)(elecProjZ)(elecEn0)
	multDistr *= maxAccelDistrFirstAp/curValAccelDistrFirstAp
endif
if(AccelDistrWaistExists)
	curValAccelDistrWaist = $nmAccelDistrWaist(elecX0)(elecZ0)
	multDistr *= maxAccelDistrWaist/curValAccelDistrWaist
endif
if(AccelDistrElecEnExists)
	curValAccelDistrElecEn = $nmAccelDistrElecEn(elecEn0)
	multDistr *= maxAccelDistrElecEn/curValAccelDistrElecEn
endif

if(absPhotEnInterv > 0)
	if(photEnIntegIntensUnit == 1)
		multDistr *= absPhotEnInterv*1000/curPhotEn
	else
		multDistr *= absPhotEnInterv*convPhW
	endif
endif

if(specTrIsDef != 0)
	multDistr *= $nmSpecTr(curPhotEn)
endif

if(multDistr != 1)
	$StoWorkName *= multDistr
endif

variable eView = srwGetSmpPhotEnStart(ObsName)
variable radCmpn2ViewAct = radCmpn2View - 1
variable xCurMin, xCurMax, yCurMin, yCurMax, rMin, rMax, r1, nxCur, nyCur
string sCurUnit

if(radCmpn2ViewAct > 0)
	SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 4, angIntens + 1, eView, xcView, zcView, 2) //vs x & y, last run
	SrwUtiGraphWindResize(10,10,210,170,0,0)
	
	xCurMin = dimoffset($nmCurWorkIntXZ, 0)
	nxCur = dimsize($nmCurWorkIntXZ, 0)
	xCurMax = xCurMin + (nxCur - 1)*dimdelta($nmCurWorkIntXZ, 0)
	sCurUnit = waveunits($nmCurWorkIntXZ, 0)
	yCurMin = dimoffset($nmCurWorkIntXZ, 1)
	nyCur = dimsize($nmCurWorkIntXZ, 1)
	yCurMax = yCurMin + (nyCur - 1)*dimdelta($nmCurWorkIntXZ, 1)
	
	if(azimAngForIntensCuts == 0)
		SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 2, angIntens + 1, eView, xcView, zcView, 2) //vs x, last run
	else //"diagonal x" cut
		rMax = (xCurMax - xcView)/cosAzimAng
		r1 = (yCurMax - zcView)/sinAzimAng
		if(rMax > r1)
			rMax = r1
		endif
		rMin = (xCurMin - xcView)/cosAzimAng
		r1 = (yCurMin - zcView)/sinAzimAng
		if(rMin < r1)
			rMin = r1
		endif
		make/O/N=(nxCur) $nmCurWorkIntX
		SetScale/I x, rMin, rMax, sCurUnit, $nmCurWorkIntX
		$nmCurWorkIntX = $nmCurWorkIntXZ(xcView + x*cosAzimAng)(zcView + x*sinAzimAng)
		display $nmCurWorkIntX; SrwUtiGraphAddFrameAndGrid()
	endif
	if(logScaleIntens)
		ModifyGraph log(left)=1
	endif
	
	SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(10,240,300,160,-1,0)
	if(absPhotEnInterv > 0)
		if(photEnIntegIntensUnit == 1)
			Label left "Ph/s/mm\\S2\\M"
		else
			Label left "W/mm\\S2\\M"
		endif
	endif

	if(azimAngForIntensCuts == 0)
		SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 3, angIntens + 1, eView, xcView, zcView, 2) //vs y, last run
	else //"diagonal y" cut
		rMin = -(xCurMin - xcView)/sinAzimAng
		r1 = (yCurMax - zcView)/cosAzimAng
		if(rMin < r1)
			rMin = r1
		endif
		rMax = -(xCurMax - xcView)/sinAzimAng
		r1 = (yCurMin - zcView)/cosAzimAng
		if(rMax > r1)
			rMax = r1
		endif
		make/O/N=(nyCur) $nmCurWorkIntZ
		SetScale/I x, rMin, rMax, sCurUnit, $nmCurWorkIntZ
		$nmCurWorkIntZ = $nmCurWorkIntXZ(xcView - x*sinAzimAng)(zcView + x*cosAzimAng)		
		display $nmCurWorkIntZ; SrwUtiGraphAddFrameAndGrid()
	endif
	if(logScaleIntens)
		ModifyGraph log(left)=1
	endif
	
	SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(10,420,300,160,-1,0)
	if(absPhotEnInterv > 0)
		if(photEnIntegIntensUnit == 1)
			Label left "Ph/s/mm\\S2\\M"
		else
			Label left "W/mm\\S2\\M"
		endif			
	endif

	SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 4, eView, xcView, zcView, 2) //vs x & y, temp. result
	SrwUtiGraphWindResize(320,10,210,170,0,0)

	xCurMin = dimoffset($nmCurStoWorkIntXZ, 0)
	nxCur = dimsize($nmCurStoWorkIntXZ, 0)
	xCurMax = xCurMin + (nxCur - 1)*dimdelta($nmCurStoWorkIntXZ, 0)
	sCurUnit = waveunits($nmCurStoWorkIntXZ, 0)
	if(angIntens)
		sCurUnit = "q"
	endif
	yCurMin = dimoffset($nmCurStoWorkIntXZ, 1)
	nyCur = dimsize($nmCurStoWorkIntXZ, 1)
	yCurMax = yCurMin + (nyCur - 1)*dimdelta($nmCurStoWorkIntXZ, 1)
	
	if(azimAngForIntensCuts == 0)
		SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 2, eView, xcView, zcView, 2) //vs x, temp. result
	else //"diagonal x" cut
		rMax = (xCurMax - xcView)/cosAzimAng
		r1 = (yCurMax - zcView)/sinAzimAng
		if(rMax > r1)
			rMax = r1
		endif
		rMin = (xCurMin - xcView)/cosAzimAng
		r1 = (yCurMin - zcView)/sinAzimAng
		if(rMin < r1)
			rMin = r1
		endif
		make/O/N=(nxCur) $nmCurStoWorkIntX
		SetScale/I x, rMin, rMax, sCurUnit, $nmCurStoWorkIntX
		$nmCurStoWorkIntX = $nmCurStoWorkIntXZ(xcView + x*cosAzimAng)(zcView + x*sinAzimAng)
		display $nmCurStoWorkIntX; SrwUtiGraphAddFrameAndGrid()
	endif
	if(logScaleIntens)
		ModifyGraph log(left)=1
	endif

	SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(320,240,300,160,-1,0)
	if(absPhotEnInterv > 0)
		if(photEnIntegIntensUnit == 1)
			Label left "Ph/s/mm\\S2\\M"
		else
			Label left "W/mm\\S2\\M"
		endif			
	endif
	
	if(azimAngForIntensCuts == 0)
		SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 3, eView, xcView, zcView, 2) //vs y, temp. result
	else //"diagonal y" cut
		rMin = -(xCurMin - xcView)/sinAzimAng
		r1 = (yCurMax - zcView)/cosAzimAng
		if(rMin < r1)
			rMin = r1
		endif
		rMax = -(xCurMax - xcView)/sinAzimAng
		r1 = (yCurMin - zcView)/cosAzimAng
		if(rMax > r1)
			rMax = r1
		endif
		make/O/N=(nyCur) $nmCurStoWorkIntZ
		SetScale/I x, rMin, rMax, sCurUnit, $nmCurStoWorkIntZ
		$nmCurStoWorkIntZ = $nmCurStoWorkIntXZ(xcView - x*sinAzimAng)(zcView + x*cosAzimAng)
		display $nmCurStoWorkIntZ; SrwUtiGraphAddFrameAndGrid()
	endif	
	if(logScaleIntens)
		ModifyGraph log(left)=1
	endif
	
	SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(320,420,300,160,-1,0)
	if(absPhotEnInterv > 0)
		if(photEnIntegIntensUnit == 1)
			Label left "Ph/s/mm\\S2\\M"
		else
			Label left "W/mm\\S2\\M"
		endif			
	endif
	DoUpdate
endif

variable multX = 0.5/(elecSigXe2*elecSigXpe2 - elecMXXp*elecMXXp)
variable BX = elecSigXe2*multX
variable GX = elecSigXpe2*multX
variable AX = elecMXXp*multX
variable SigPX = 1/sqrt(2*GX)
variable SigQX = sqrt(GX/(2*(BX*GX - AX*AX)))

variable multZ = 0.5/(elecSigZe2*elecSigZpe2 - elecMZZp*elecMZZp)
variable BZ = elecSigZe2*multZ
variable GZ = elecSigZpe2*multZ
variable AZ = elecMZZp*multZ
variable SigPZ = 1/sqrt(2*GZ)
variable SigQZ = sqrt(GZ/(2*(BZ*GZ - AZ*AZ)))

variable auxPX, auxPXp, auxPZ, auxPZp
variable xStartWfr, xEndWfr, xNpWfr, xStepWfr, zStartWfr, zEndWfr, zNpWfr, zStepWfr

make/O $nmXcSigLPTau = {{0,0,0,0,0},{1,1,1,1,1}}
make/O/N=5 $nmElecInitCond
variable/C auxPosAng
variable elecX, elecXp, elecZ, elecZp, elecEn
variable xbeg, xstep, xend, nx, zbeg, zstep, zend, nz
variable nxExtraLeft, nxExtraRight, nxExtra, nzExtraLeft, nzExtraRight, nzExtra
variable iElec = 2, inv_iElec
variable perSave = 10, iSave = 0, initRand = 1, elecCoordAreInLimits
variable modeRand = 0 // 0- random, 1- LPTau

do
	if(treatElecEnergySpreadOnly)
		elecEn = elecEn0 + gnoise(elecAbsEnSpr)
		elecX = elecX
		elecXp = elecXp0	
		elecZ = elecZ0
		elecZp = elecZp0
	else
		$nmElecInitCond = 0
		srUtiRandGsn($nmXcSigLPTau, 5, initRand, modeRand, $nmElecInitCond)
		initRand = 0
	
		auxPXp = SigQX*$nmElecInitCond[0]
		auxPX = SigPX*$nmElecInitCond[1] + AX*auxPXp/GX
		elecX = elecX0 + auxPX
		elecXp = elecXp0 + auxPXp
	
		auxPZp = SigQZ*$nmElecInitCond[2]
		auxPZ = SigPZ*$nmElecInitCond[3] + AZ*auxPZp/GZ
		elecZ = elecZ0 + auxPZ
		elecZp = elecZp0 + auxPZp
	
		elecEn = elecEn0 + elecAbsEnSpr*$nmElecInitCond[4]
	endif
	
	//$nmElecInitCond = 0
	//srUtiRandGsn($nmXcSigLPTau, 5, initRand, modeRand, $nmElecInitCond)
	//initRand = 0
	//
	//auxPXp = SigQX*$nmElecInitCond[0]
	//auxPX = SigPX*$nmElecInitCond[1] + AX*auxPXp/GX
	//elecX = elecX0 + auxPX
	//elecXp = elecXp0 + auxPXp
	//
	//auxPZp = SigQZ*$nmElecInitCond[2]
	//auxPZ = SigPZ*$nmElecInitCond[3] + AZ*auxPZp/GZ
	//elecZ = elecZ0 + auxPZ
	//elecZp = elecZp0 + auxPZp
	//
	//elecEn = elecEn0 + elecAbsEnSpr*$nmElecInitCond[4]
	
	$StoWorkName *= (iElec - 1)/iElec

	elecCoordAreInLimits = 1
	if((elecEn < elecLimitEstart) %| (elecEn > elecLimitEend))
		elecCoordAreInLimits = 0
	endif
	if((elecCoordAreInLimits != 0) %& ((elecX < elecLimitXstart) %| (elecX > elecLimitXend)))
		elecCoordAreInLimits = 0
	endif
	if((elecCoordAreInLimits != 0) %& ((elecXp < elecLimitXPstart) %| (elecXp > elecLimitXPend)))
		elecCoordAreInLimits = 0
	endif
	if((elecCoordAreInLimits != 0) %& ((elecZ < elecLimitZstart) %| (elecZ > elecLimitZend)))
		elecCoordAreInLimits = 0
	endif
	if((elecCoordAreInLimits != 0) %& ((elecZp < elecLimitZPstart) %| (elecZp > elecLimitZPend)))
		elecCoordAreInLimits = 0
	endif
	propagWasDone = 0
	AccelDistrFirstApPasses = 1
	AccelDistrWaistPasses = 1
	AccelDistrElecEnPasses = 1
	multDistr = 1
	
	if(elecCoordAreInLimits != 0)
		
		if(AccelDistrFirstApExists)
			elecProjX = -elecX - rProj*elecXp
			elecProjZ = -elecZ - rProj*elecZp
			
			curValAccelDistrFirstAp = -1
			if((elecEn > minElecEnAccelDistrFirstAp) %& (elecEn < maxElecEnAccelDistrFirstAp))
				if((elecProjX > minProjXAccelDistrFirstAp) %& (elecProjX < maxProjXAccelDistrFirstAp))
					if((elecProjZ > minProjZAccelDistrFirstAp) %& (elecProjZ < maxProjZAccelDistrFirstAp))
						curValAccelDistrFirstAp = Interp3d($nmAccelDistrFirstAp, elecProjX, elecProjZ, elecEn)
					endif
				endif
			endif
			if(curValAccelDistrFirstAp < 0)
				curValAccelDistrFirstAp = $nmAccelDistrFirstAp(elecProjX)(elecProjZ)(elecEn)
			endif
			testValAccelDistrFirstAp = halfMaxAccelDistrFirstAp + enoise(halfMaxAccelDistrFirstAp)
			if((curValAccelDistrFirstAp <= 0) %| (curValAccelDistrFirstAp <= testValAccelDistrFirstAp))
				AccelDistrFirstApPasses = 0
			else
				multDistr *= maxAccelDistrFirstAp/curValAccelDistrFirstAp
			endif
		endif
		
		if(AccelDistrElecEnExists)
			testValAccelDistrElecEn = halfMaxAccelDistrElecEn + enoise(halfMaxAccelDistrElecEn)
			curValAccelDistrElecEn = $nmAccelDistrElecEn(elecEn)
			if((curValAccelDistrElecEn <= 0) %| (curValAccelDistrElecEn < testValAccelDistrElecEn))
				AccelDistrElecEnPasses = 0
			else
				multDistr *= maxAccelDistrElecEn/curValAccelDistrElecEn
			endif
		endif
		
		if(AccelDistrFirstApPasses %& AccelDistrWaistExists)
			curValAccelDistrWaist = -1
			if((elecX > minXAccelDistrWaist) %& (elecX < maxXAccelDistrWaist))
				if((elecZ > minZAccelDistrWaist) %& (elecZ < maxZAccelDistrWaist))
					curValAccelDistrWaist = Interp2d($nmAccelDistrWaist, elecX, elecZ)
				endif
			endif
			if(curValAccelDistrWaist < 0)
				curValAccelDistrWaist = $nmAccelDistrWaist(elecX)(elecZ)
			endif
			testValAccelDistrWaist = halfMaxAccelDistrWaist + enoise(halfMaxAccelDistrWaist)
			if((curValAccelDistrWaist <= 0) %| (curValAccelDistrWaist <= testValAccelDistrWaist))
				AccelDistrWaistPasses = 0
			else
				multDistr *= maxAccelDistrWaist/curValAccelDistrWaist
			endif
		endif
		
		if(AccelDistrFirstApPasses %& AccelDistrWaistPasses %& AccelDistrElecEnPasses)
	
			print "i=", iElec, "  Electron Coordinates: x=", elecX, " x'=", elecXp, "z=", elecZ, " z'=", elecZp, "energy=", elecEn
	
			SrwElecFilament(ElecWorkCore, elecEn, elecCur, elecS0, elecX*1000, elecXp*1000, elecZ*1000, elecZp*1000)
			
			SrwUtiTriggerPrint(2)
			SrwElecThick(ElecWorkName,0,0,0,0,0,0,0,0,0)
			SrwUtiTriggerPrint(1)
			
			if(absPhotEnInterv > 0)
				curPhotEn = cenPhotEn + enoise(halfAbsPhotEnInterv)
				srwSetSmpPhotEnStart(ObsName, curPhotEn)
				srwSetSmpPhotEnEnd(ObsName, curPhotEn)
			endif
			SrwWfrCreate(RadWorkCore, ElecWorkName, MagName, ObsName, ObsNxNzForProp, abs(ObsNxNzSamplFact))
	
			//estimate flux, and if it's too small - don't propagate
			SrwWfr2Int(RadWorkName,sufRadWorkFlux,7,5,1,1,1,0,0,1)
			curFlux = $nmRadWorkFlux[0]
	
			if(curFlux >= relFluxRejectPropag*curMaxFlux)
				
				SrwWfrPropList(RadWorkName, PropBL, 1, "")
			
				//estimate propagated flux: to check for a possible errror
				//SrwWfr2Int(RadWorkName,sufRadWorkFlux,7,5,1,1,1,0,0,1)
				//curPropFlux = $nmRadWorkFlux[0]
				
				if(angIntens)
					srWfrSetRepres($RadWorkName, "A")
				endif
				
				xStartWfr = dimoffset($nmRadWorkEX, 1)
				xNpWfr = dimsize($nmRadWorkEX, 1)	
				xStepWfr = dimdelta($nmRadWorkEX, 1)
				xEndWfr = xStartWfr + (xNpWfr - 1)*xStepWfr
				zStartWfr = dimoffset($nmRadWorkEX, 2)
				zNpWfr = dimsize($nmRadWorkEX, 2)	
				zStepWfr = dimdelta($nmRadWorkEX, 2)
				zEndWfr = zStartWfr + (zNpWfr - 1)*zStepWfr
				
				multDistr /= iElec
				
				if(absPhotEnInterv != 0)
					if(photEnIntegIntensUnit == 1)
						multDistr *= absPhotEnInterv*1000/curPhotEn
					else
						multDistr *= absPhotEnInterv*convPhW
					endif
				endif
				
				if(specTrIsDef != 0)
					multDistr *= $nmSpecTr(curPhotEn)
				endif
				
				$StoWorkName += multDistr*srwUtiE2Stokes($nmRadWorkEX(y)(z)(t), $nmRadWorkEZ(y)(z)(t), p)*srwUtiNonZeroIntervB(z, xStartWfr, xEndWfr)*srwUtiNonZeroIntervB(t, zStartWfr, zEndWfr) //assuming both components exist
				propagWasDone = 1
				
			else
				propagWasDone = 0
			endif
			if(curMaxFlux < curFlux)
				curMaxFlux = curFlux
			endif

			if((radCmpn2ViewAct > 0) %& (propagWasDone > 0))
				SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 4, angIntens + 1, eView, xcView, zcView, 1) //vs x & y, last run
				
				xCurMin = dimoffset($nmCurWorkIntXZ, 0)
				nxCur = dimsize($nmCurWorkIntXZ, 0)
				xCurMax = xCurMin + (nxCur - 1)*dimdelta($nmCurWorkIntXZ, 0)
				sCurUnit = waveunits($nmCurWorkIntXZ, 0)
				yCurMin = dimoffset($nmCurWorkIntXZ, 1)
				nyCur = dimsize($nmCurWorkIntXZ, 1)
				yCurMax = yCurMin + (nyCur - 1)*dimdelta($nmCurWorkIntXZ, 1)
				
				if(azimAngForIntensCuts == 0)
					SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 2, angIntens + 1, eView, xcView, zcView, 1) //vs x, last run
					SrwWfr2Int(RadWorkName, IntCurSuf, radCmpn2ViewAct, 1, 3, angIntens + 1, eView, xcView, zcView, 1) //vs y, last run
				else
					rMax = (xCurMax - xcView)/cosAzimAng
					r1 = (yCurMax - zcView)/sinAzimAng
					if(rMax > r1)
						rMax = r1
					endif
					rMin = (xCurMin - xcView)/cosAzimAng
					r1 = (yCurMin - zcView)/sinAzimAng
					if(rMin < r1)
						rMin = r1
					endif
					make/O/N=(nxCur) $nmCurWorkIntX
					SetScale/I x, rMin, rMax, sCurUnit, $nmCurWorkIntX
					$nmCurWorkIntX = $nmCurWorkIntXZ(xcView + x*cosAzimAng)(zcView + x*sinAzimAng)

					rMin = -(xCurMin - xcView)/sinAzimAng
					r1 = (yCurMax - zcView)/cosAzimAng
					if(rMin < r1)
						rMin = r1
					endif
					rMax = -(xCurMax - xcView)/sinAzimAng
					r1 = (yCurMin - zcView)/cosAzimAng
					if(rMax > r1)
						rMax = r1
					endif
					make/O/N=(nyCur) $nmCurWorkIntZ
					SetScale/I x, rMin, rMax, sCurUnit, $nmCurWorkIntZ
					$nmCurWorkIntZ = $nmCurWorkIntXZ(xcView - x*sinAzimAng)(zcView + x*cosAzimAng)
				endif

				SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 4, eView, xcView, zcView, 1) //vs x & y, temp. result
				
				xCurMin = dimoffset($nmCurStoWorkIntXZ, 0)
				nxCur = dimsize($nmCurStoWorkIntXZ, 0)
				xCurMax = xCurMin + (nxCur - 1)*dimdelta($nmCurStoWorkIntXZ, 0)
				sCurUnit = waveunits($nmCurStoWorkIntXZ, 0)
				if(angIntens)
					sCurUnit = "q"
				endif
				yCurMin = dimoffset($nmCurStoWorkIntXZ, 1)
				nyCur = dimsize($nmCurStoWorkIntXZ, 1)
				yCurMax = yCurMin + (nyCur - 1)*dimdelta($nmCurStoWorkIntXZ, 1)
				
				if(azimAngForIntensCuts == 0)
					SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 2, eView, xcView, zcView, 1) //vs x, temp. result
					SrwSto2Int(StoWorkName, IntResSuf, radCmpn2ViewAct, 3, eView, xcView, zcView, 1) //vs y, temp. result
				else
					rMax = (xCurMax - xcView)/cosAzimAng
					r1 = (yCurMax - zcView)/sinAzimAng
					if(rMax > r1)
						rMax = r1
					endif
					rMin = (xCurMin - xcView)/cosAzimAng
					r1 = (yCurMin - zcView)/sinAzimAng
					if(rMin < r1)
						rMin = r1
					endif
					make/O/N=(nxCur) $nmCurStoWorkIntX
					SetScale/I x, rMin, rMax, sCurUnit, $nmCurStoWorkIntX
					$nmCurStoWorkIntX = $nmCurStoWorkIntXZ(xcView + x*cosAzimAng)(zcView + x*sinAzimAng)

					rMin = -(xCurMin - xcView)/sinAzimAng
					r1 = (yCurMax - zcView)/cosAzimAng
					if(rMin < r1)
						rMin = r1
					endif
					rMax = -(xCurMax - xcView)/sinAzimAng
					r1 = (yCurMin - zcView)/cosAzimAng
					if(rMax > r1)
						rMax = r1
					endif
					make/O/N=(nyCur) $nmCurStoWorkIntZ
					SetScale/I x, rMin, rMax, sCurUnit, $nmCurStoWorkIntZ
					$nmCurStoWorkIntZ = $nmCurStoWorkIntXZ(xcView - x*sinAzimAng)(zcView + x*cosAzimAng)
				endif
				DoUpdate
			endif
		endif
	endif
	
	iElec += 1
	if(propagWasDone)
		iSave += 1
	endif	
	if(iSave >= perSave)
		SrwWfrDel(RadWorkName)
		saveexperiment
		iSave = 0
	endif
	
while(iElec <= MaxPrt)

KillWaves/Z  $ElecWorkName, $nmXcSigLPTau, $nmElecInitCond
end

//+++++++++++++++++++++++++++++++++++++++
//Auxiliary input proc for multi-e wavefront propagation
//(SrwWfrEmitPropStokesMultiE)
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrEmitPropStokesAuxIn(initWfrCalcMeth, initWfrCalcPrec, angIntens, azimAngForIntensCuts, scaleIntens, relFluxRejectProp, treatElecEnSpreadOnly, photEnIntegIntensUnit, nmSpecTr)
variable angIntens = srwUtiGetValN("angIntens", 1, "SrwWfrEmitPropStokesAuxIn") 
variable azimAngForIntensCuts = srwUtiGetValN("azimAngForIntensCuts", 0, "SrwWfrEmitPropStokesAuxIn") 
variable scaleIntens = srwUtiGetValN("scaleIntens", 1, "SrwWfrEmitPropStokesAuxIn")
variable initWfrCalcMeth = SrwMode
variable initWfrCalcPrec = srwUtiGetValN("initWfrCalcPrec", SrwPrec, "SrwWfrEmitPropStokesAuxIn")
variable relFluxRejectProp = srwUtiGetValN("relFluxRejectProp", 1e-04, "SrwWfrEmitPropStokesAuxIn")
variable treatElecEnSpreadOnly = srwUtiGetValN("treatElecEnSpreadOnly", 1, "SrwWfrEmitPropStokesAuxIn")
variable photEnIntegIntensUnit = srwUtiGetValN("photEnIntegIntensUnit", 1, "SrwWfrEmitPropStokesAuxIn")
string nmSpecTr = srwUtiGetValS("nmSpecTr", "", "SrwWfrEmitPropStokesAuxIn")
prompt angIntens, "Coordinate or Angular Final Intens. Repres.", popup "Coordinate;Angular"
prompt azimAngForIntensCuts, "Azimuthal Angle for Intensity Cuts [rad]"
prompt scaleIntens, "Extracted Intensity Scale", popup "Linear;Logarithmic"
prompt initWfrCalcMeth, "Initial Wavefront Calculation Method", popup SrwPOPUPMode
prompt initWfrCalcPrec, "Integration Step or Relative Precision"
prompt relFluxRejectProp, "Rel. Flux Threshold to Skip Propagation"
prompt treatElecEnSpreadOnly, "Treat Energy Spread Only, No Emittance", popup "No;Yes"
prompt photEnIntegIntensUnit, "Intensity Units (if integ. vs Phot. En.)", popup "Ph/s/mm^2;W/mm^2"
prompt nmSpecTr, "Spectral Transmission wave", popup  "_none_;" + Wavelist("*",";","TEXT:0,DIMS:1")
Silent 1						|	Propagating Wavefront ...
PauseUpdate

SrwMode = initWfrCalcMeth
if(initWfrCalcMeth == 1)
	SrwRadIntStep = initWfrCalcPrec
else
	SrwPrec = initWfrCalcPrec
endif
srwUtiSetValN("angIntens", angIntens, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("azimAngForIntensCuts", azimAngForIntensCuts, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("scaleIntens", scaleIntens, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("initWfrCalcPrec", initWfrCalcPrec, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("relFluxRejectProp", relFluxRejectProp, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("treatElecEnSpreadOnly", treatElecEnSpreadOnly, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("photEnIntegIntensUnit", photEnIntegIntensUnit, "SrwWfrEmitPropStokesAuxIn")
srwUtiSetValN("justCalled", 1, "SrwWfrEmitPropStokesAuxIn")
if(cmpstr(nmSpecTr, "_none_") == 0)
	nmSpecTr = ""
endif
srwUtiSetValS("nmSpecTr", nmSpecTr, "SrwWfrEmitPropStokesAuxIn")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Propagate wavefront to deduce Stokes params of radiation emitted by
//Thick electron beam.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropagStokesMultiE(StoName, ElecName, Wfr, BL, PrecParM, MaxPrt)
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

//+++++++++++++++++++++++++++++++++++++++
//
//Propagate wavefront to deduce Stokes params of radiation emitted by
//Thick electron beam.
//In this version, initial wavefront is re-calculated for every macro-particle.
//Old version
//+++++++++++++++++++++++++++++++++++++++
//proc SrwWfrEmitPropagStokesMultiE(StoName, ElecName, MagName, ObsName, ObsNxNzForProp, ObsNxNzSamplFactX, ObsNxNzSamplFactZ, BL, MaxPrt)
proc SrwWfrEmitPropagStokesMultiE(StoName, ElecName, MagName, ObsName, ObsNxNzForProp, ObsNxNzSamplFact, BL, MaxPrt, xcView, zcView)
string StoName=srwUtiGetValS("StoName", "Stk", "SrwWfrEmitPropagStokesMultiE")
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagName+SrwFieldType
string ObsName=SrwSmpName+SrwSmpType
variable ObsNxNzForProp=SrwSmpNxNzForProp
//variable ObsNxNzSamplFactX=SrwSmpNxNzSamplFact
//variable ObsNxNzSamplFactZ=SrwSmpNxNzSamplFact
variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
string BL=SrwBliLast+SrwBeamlineType
//variable PrecParM=srwUtiGetValN("PrecParM", 1, "SrwWfrPropagStokesMultiE")
variable MaxPrt=srwUtiGetValN("MaxPrt", 10000, "SrwWfrEmitPropagStokesMultiE") 
variable xcView=srwUtiGetValN("xcView", 0, "SrwWfrEmitPropagStokesMultiE") 
variable zcView=srwUtiGetValN("zcView", 0, "SrwWfrEmitPropagStokesMultiE") 
prompt StoName, "Name of the Stokes structure"
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
//prompt ObsNxNzSamplFactX,"Oversampling Factor in Hor. Direction"
//prompt ObsNxNzSamplFactZ,"Oversampling Factor in Vert. Direction" //SrwPSmpNxNzSamplFact
prompt ObsNxNzSamplFact,"Oversampling Factor" //SrwPSmpNxNzSamplFact
prompt BL, "Optical Element",popup Wavelist("*"+SrwBeamlineType ,";", "")
//prompt PrecParM, "Multi-e propag. precision"
prompt xcView, "Hor. Center Point for Viewing [mm]"
prompt zcView, "Vert. Center Point for Viewing [mm]"
Silent 1						|	Propagating the Wavefront ...
//PauseUpdate

SrwBliLast = BL[0,strlen(BL)-strlen(SrwBeamlineType)-1]
SrwElecName = ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName = ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwSmpNxNzForProp = ObsNxNzForProp
SrwSmpNxNzSamplFact = ObsNxNzSamplFact

//SrwRadName = Wfr[0,strlen(Wfr)-strlen(SrwRadType)-1]
srwUtiSetValS("StoName", StoName, "SrwWfrEmitPropagStokesMultiE")
//srwUtiSetValN("PrecParM", PrecParM, "SrwWfrPropagStokesMultiE")
srwUtiSetValN("MaxPrt", MaxPrt, "SrwWfrEmitPropagStokesMultiE")
srwUtiSetValN("xcView", xcView, "SrwWfrEmitPropagStokesMultiE")
srwUtiSetValN("zcView", zcView, "SrwWfrEmitPropagStokesMultiE")

variable elecSigX = srwGetElecBeamHorSizeRMS(ElecName)
variable elecSigXp = srwGetElecBeamHorDivergRMS(ElecName)
variable elecMXXp = srwGetElecBeamHorMixedMom(ElecName)
variable elecSigXe2 = elecSigX*elecSigX
variable elecSigXpe2 = elecSigXp*elecSigXp
variable elecSigZ = srwGetElecBeamVertSizeRMS(ElecName)
variable elecSigZp = srwGetElecBeamVertDivergRMS(ElecName)
variable elecMZZp = srwGetElecBeamVertMixedMom(ElecName)
variable elecSigZe2 = elecSigZ*elecSigZ
variable elecSigZpe2 = elecSigZp*elecSigZp
variable elecRelEnSpr = srwGetElecBeamRelEnSprRMS(ElecName)

variable elecEn0 = srwGetElecBeamEnergy(ElecName)
variable elecAbsEnSpr = elecEn0*elecRelEnSpr

variable elecCur = srwGetElecBeamCurrent(ElecName)
variable elecS0 = srwGetElecBeamLongPos(ElecName)
variable elecX0 = srwGetElecBeamHorPos(ElecName)
variable elecXp0 = srwGetElecBeamHorAng(ElecName)
variable elecZ0 = srwGetElecBeamVertPos(ElecName)
variable elecZp0 = srwGetElecBeamVertAng(ElecName)

string ElecWorkCore = "elecEmitPropStokesMultiE"
SrwElecFilament(ElecWorkCore, ElecEn0, ElecCur, elecS0, elecX0*1000, elecXp0*1000, elecZ0*1000, elecZp0*1000)
string ElecWorkName = ElecWorkCore + "_ebm"
string RadWorkCore = "wfrEmitPropStokesMultiE"
SrwWfrCreate(RadWorkCore, ElecWorkName, MagName, ObsName, ObsNxNzForProp, ObsNxNzSamplFact)
string RadWorkName = RadWorkCore + "_rad"
variable UseResBefore = 1, UseResAfter = 1, PrecParamProp = 1, DplBeforeProp = 1
string RadDplName = "dummy"

	//Hack:
	SrwWfrResize(RadWorkName,1,1,1,1,10,1,"Wfrd")
	//SrwWfrResize(RadWorkName,1,1,50,1,10,1,"Wfrd")
	//End hack

SrwWfrPropag(RadWorkName, BL, UseResBefore, UseResAfter, PrecParamProp, DplBeforeProp, RadDplName)

SrwWfr2Int(RadWorkName, "Work", 7, 1, 4, 1, 1, 0, 0, 2)
SrwUtiGraphWindResize(320,10,250,220,0,0)

string IntName = RadWorkCore + "Work_xz"

duplicate/O $IntName $StoName
variable xbeg0 = dimoffset($StoName, 0)
variable xstep0 = dimdelta($StoName, 0)
variable nx0 = dimsize($StoName, 0)
variable xend0 = xbeg0 + (nx0 - 1)*xstep0
variable zbeg0 = dimoffset($StoName, 1)
variable zstep0 = dimdelta($StoName, 1)
variable nz0 = dimsize($StoName, 1)
variable zend0 = zbeg0 + (nz0 - 1)*zstep0

string StoNameCutVsX = StoName + "_cutX", StoNameCutVsZ = StoName + "_cutZ"
make/O/N=(nx0) $StoNameCutVsX
SetScale/P x xbeg0,xstep0,"", $StoNameCutVsX
make/O/N=(nz0) $StoNameCutVsZ
SetScale/P x zbeg0,zstep0,"", $StoNameCutVsZ
//variable xcView = 0.5*(xbeg0 + xend0), zcView = 0.5*(zbeg0 + zend0)
$StoNameCutVsX = $StoName(x)(zcView*0.001)
$StoNameCutVsZ = $StoName(xcView*0.001)(x)

variable xbeg0_old = xbeg0, xend0_old = xend0
variable zbeg0_old = zbeg0, zend0_old = zend0
string AuxStoName = "auxEmitPropStokesMultiE"

Display; AppendImage $StoName; SrwUtiGraphWindResize(10,10,250,220,0,0)
Display $StoNameCutVsX; SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(10,320,300,200,-1,0)
Display $StoNameCutVsZ; SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(320,320,300,200,-1,0)
DoUpdate

variable/C auxPosAng
variable elecX, elecXp, elecZ, elecZp, elecEn
variable xbeg, xstep, xend, nx, zbeg, zstep, zend, nz
variable nxExtraLeft, nxExtraRight, nxExtra, nzExtraLeft, nzExtraRight, nzExtra
variable iElec = 2, inv_iElec
variable perSave = 10, iSave = 0
do
	auxPosAng = srwUtiRandGauss2D(elecSigXe2, elecMXXp, elecSigXpe2)
	elecX = elecX0 + real(auxPosAng)
	elecXp = elecXp0 + imag(auxPosAng)
	auxPosAng = srwUtiRandGauss2D(elecSigZe2, elecMZZp, elecSigZpe2)
	elecZ = elecZ0 + real(auxPosAng)
	elecZp = elecZp0 + imag(auxPosAng)
	elecEn = elecEn0 + gnoise(elecAbsEnSpr)
	
		print "x=", elecX, " xp=", elecXp, "z=", elecZ, " zp=", elecZp, "en=", elecEn
	
	SrwElecFilament(ElecWorkCore, elecEn, ElecCur, elecS0, elecX*1000, elecXp*1000, elecZ*1000, elecZp*1000)
	SrwUtiTriggerPrint(2)
	SrwElecThick(ElecWorkCore + "_ebm",0,0,0,0,0,0,0,0,0)
	SrwUtiTriggerPrint(1)
	SrwWfrCreate(RadWorkCore, ElecWorkName, MagName, ObsName, ObsNxNzForProp, ObsNxNzSamplFact)
	
		//Hack:
		SrwWfrResize(RadWorkName,1,1,1,1,10,1,"Wfrd")
		//SrwWfrResize(RadWorkName,1,1,50,1,10,1,"Wfrd")
		//End hack

	SrwWfrPropag(RadWorkName, BL, UseResBefore, UseResAfter, PrecParamProp, DplBeforeProp, RadDplName)
	
	SrwWfr2Int(RadWorkName, "Work", 7, 1, 4, 1, 1, 0, 0, 1)
	
	xbeg = dimoffset($IntName, 0)
	xstep = dimdelta($IntName, 0)
	nx = dimsize($IntName, 0)
	xend = xbeg + (nx - 1)*xstep
	nxExtraLeft = 0; nxExtraRight = 0
	if(xbeg < xbeg0)
		nxExtraLeft = trunc((xbeg0 - xbeg)/xstep0)
	endif
	if(xend > xend0)
		nxExtraRight = trunc((xend - xend0)/xstep0)
	endif
	nxExtra = nxExtraLeft + nxExtraRight
	
	zbeg = dimoffset($IntName, 1)
	zstep = dimdelta($IntName, 1)
	nz = dimsize($IntName, 1)
	zend = zbeg + (nz - 1)*zstep
	nzExtraLeft = 0; nzExtraRight = 0
	if(zbeg < zbeg0)
		nzExtraLeft = trunc((zbeg0 - zbeg)/zstep0)
	endif
	if(zend > zend0)
		nzExtraRight = trunc((zend - zend0)/zstep0)
	endif
	nzExtra = nzExtraLeft + nzExtraRight

	if((nxExtra > 0) %| (nzExtra > 0))
		//if((nxExtra > 0.25*nx) %| (nzExtra > 0.25*nz))
		//	print "Requirement of abnormal wavefront increase ignored"
		//	//continue
		//else
			duplicate/O $StoName $AuxStoName
			if(nxExtra > 0)
				nx0 += nxExtra
				if(nxExtraLeft > 0)
					xbeg0_old = xbeg0
					xbeg0 -= nxExtraLeft*xstep0
				endif
				xend0_old = xend0
				xend0 = xbeg0 + (nx0 - 1)*xstep0
			endif
			if(nzExtra > 0)
				nz0 += nzExtra
				if(nzExtraLeft > 0)
					zbeg0_old = zbeg0
					zbeg0 -= nzExtraLeft*zstep0
				endif
				zend0_old = zend0
				zend0 = zbeg0 + (nz0 - 1)*zstep0
			endif
		
			redimension/N=(nx0, nz0) $StoName
			SetScale/P x xbeg0,xstep0,"", $StoName
			SetScale/P y zbeg0,zstep0,"", $StoName
			$StoName = 0
			$StoName += $AuxStoName(x)(y)*srwUtiNonZeroInterval(x, xbeg0_old, xend0_old)*srwUtiNonZeroInterval(y, zbeg0_old, zend0_old)
		
			redimension/N=(nx0) $StoNameCutVsX
			SetScale/P x xbeg0,xstep0,"", $StoNameCutVsX
			redimension/N=(nz0) $StoNameCutVsZ
			SetScale/P x zbeg0,zstep0,"", $StoNameCutVsZ
		//endif
	endif
	
	//if((nxExtra <= 0.5*nx) %& (nzExtra <= 0.5*nz))
		inv_iElec = 1/iElec
		$StoName *= (iElec - 1)*inv_iElec
		$StoName += inv_iElec*$IntName(x)(y)*srwUtiNonZeroInterval(x, xbeg, xend)*srwUtiNonZeroInterval(y, zbeg, zend)
	
		//xcView = 0.5*(xbeg0 + xend0); zcView = 0.5*(zbeg0 + zend0)
		$StoNameCutVsX = $StoName(x)(zcView*0.001)
		$StoNameCutVsZ = $StoName(xcView*0.001)(x)
		DoUpdate

		iElec += 1
	//endif
	
	iSave += 1
	if(iSave >= perSave)
		saveexperiment
		iSave = 0
	endif
	
while(iElec <= MaxPrt)

KillWaves/Z  $ElecWorkName, $IntName
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculate SR Intensity Distributions for different Electron Energies
//within energy spread of a Thick Electron Beam.
//Auxiliary function; to be used for acceleration of partially-coherent 
//wavefronts propagation simulations.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrEmitPropIntDistrib(IntName, ElecName, MagName, ObsName, RadCmpn, ElecEnRangeSig, NumElecEn) //, Dx, Dz)
string IntName=srwUtiGetValS("IntName", "Int", "SrwWfrEmitPropIntDistrib")
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagName+SrwFieldType
string ObsName=SrwSmpName+SrwSmpType
variable RadCmpn=srwUtiGetValN("RadCmpn", 7, "SrwWfrEmitPropIntDistrib")
variable ElecEnRangeSig=srwUtiGetValN("ElecEnRangeSig", 4, "SrwWfrEmitPropIntDistrib")
variable NumElecEn=srwUtiGetValN("NumElecEn", 11, "SrwWfrEmitPropIntDistrib")
//variable Dx=srwUtiGetValN("Dx", 1, "SrwWfrEmitPropIntDistrib")
//variable Dz=srwUtiGetValN("Dz", 1, "SrwWfrEmitPropIntDistrib")
prompt IntName, "Name for the Resulting Distribution"
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt RadCmpn, "Polarization Component to consider", popup SrwPOPUPPolar+";Total"
prompt ElecEnRangeSig, "Electron Energy Range in Stand. Dev."
prompt NumElecEn, "Number of Electron Energies"
//prompt Dx, "Horizontal Aperture Size [mm]"
//prompt Dz, "Vertical Aperture Size [mm]"
Silent 1						|	Calculating SR ...
PauseUpdate

SrwElecName = ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName = ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]

srwUtiSetValS("IntName", IntName, "SrwWfrEmitPropIntDistrib")
srwUtiSetValN("RadCmpn", RadCmpn, "SrwWfrEmitPropIntDistrib")
srwUtiSetValN("ElecEnRangeSig", ElecEnRangeSig, "SrwWfrEmitPropIntDistrib")
srwUtiSetValN("NumElecEn", NumElecEn, "SrwWfrEmitPropIntDistrib")

//srwUtiSetValN("Dx", Dx, "SrwWfrEmitPropIntDistrib")
//srwUtiSetValN("Dz", Dz, "SrwWfrEmitPropIntDistrib")

string ElecWorkCore = "elecEmitPropIntDistrib"
string ElecWorkName = ElecWorkCore + "_ebm"
string RadWorkCore = "wfrEmitPropIntDistrib"
string RadWorkName = RadWorkCore + "_rad"
string IntWorkSuf = "Icur"
string IntWorkName = RadWorkCore + IntWorkSuf + "_xz"

variable radSampNpX = srwGetSmpHorPosNp(ObsName)
variable radSampStartX = srwGetSmpHorPosStart(ObsName)
variable radSampEndX = srwGetSmpHorPosEnd(ObsName)
variable radSampNpZ = srwGetSmpVertPosNp(ObsName)
variable radSampStartZ = srwGetSmpVertPosStart(ObsName)
variable radSampEndZ = srwGetSmpVertPosEnd(ObsName)
//variable radSampStepX = 0,  radSampStepZ = 0
//if(radSampNpX > 1)
//	radSampStepX = (radSampEndX - radSampStartX)/(radSampNpX - 1)
//endif
//if(radSampNpZ > 1)
//	radSampStepZ = (radSampEndZ - radSampStartZ)/(radSampNpZ - 1)
//endif

variable photEn_keV = 0.001*srwGetSmpPhotEnStart(ObsName)

variable elecCur = srwGetElecBeamCurrent(ElecName)
variable elecS0 = srwGetElecBeamLongPos(ElecName)
variable elecX0 = srwGetElecBeamHorPos(ElecName)
variable elecXp0 = srwGetElecBeamHorAng(ElecName)
variable elecZ0 = srwGetElecBeamVertPos(ElecName)
variable elecZp0 = srwGetElecBeamVertAng(ElecName)

variable elecEn0 = srwGetElecBeamEnergy(ElecName)
variable elecEnRelSpr = srwGetElecBeamRelEnSprRMS(ElecName)
variable elecEnAbsSpr = elecEn0*elecEnRelSpr
variable elecEnRange = ElecEnRangeSig*elecEnAbsSpr
variable elecEnStep = 0
if(NumElecEn > 1)
	elecEnStep = elecEnRange/(NumElecEn - 1)
endif

variable elecEnStart = elecEn0 - 0.5*elecEnRange
variable elecEnEnd = elecEn0 + 0.5*elecEnRange

if(NumElecEn > 1)
	make/O/N=(radSampNpX, radSampNpZ, NumElecEn) $IntName
	SetScale/I z elecEnStart, elecEnEnd, "GeV", $IntName
else
	make/O/N=(radSampNpX, radSampNpZ) $IntName
endif
SetScale/I x radSampStartX, radSampEndX, "m", $IntName
SetScale/I y radSampStartZ, radSampEndZ, "m", $IntName

variable elecEn = elecEn0 - 0.5*elecEnRange
variable iElec = 0

SrwUtiTriggerPrint(2)

do
	print "Electron Energy: ", elecEn

	SrwElecFilament(ElecWorkCore, elecEn, elecCur, elecS0, elecX0*1000, elecXp0*1000, elecZ0*1000, elecZp0*1000)
	SrwElecThick(ElecWorkName,0,0,0,0,0,0,0,0,0)
	SrwWfrCreate(RadWorkCore, ElecWorkName, MagName, ObsName, 1, 1)
	SrwWfr2Int(RadWorkName, IntWorkSuf, RadCmpn, 1, 4, 1, photEn_keV, 0, 0, 1)
	
	if(NumElecEn > 1)
		$IntName[][][iElec] = $IntWorkName[p][q]
	else
		$IntName = $IntWorkName[p][q]
	endif

	elecEn += elecEnStep
	iElec += 1
while(iElec < NumElecEn)

SrwUtiTriggerPrint(1)
//killwaves/Z $IntWorkName
//SrwWfrDel(RadWorkName)
end