
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculate Intensity of radiation taking into account electron beam emittance and energy spread
//by summing-up contributions from different "macro-electrons"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwIntArbMCCreate(RadName, ElecName, MagName, ObsName, RadCmpnType, NumPt, Mode1D, Prec1D, MethTransvEmit)
string RadName=srwUtiGetValS("RadName", "IntME", "SrwIntArbMCCreate")
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string ObsName=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
variable RadCmpnType=SrwViewRadCmpnType
variable NumPt=srwUtiGetValN("NumPt", 100, "SrwIntArbMCCreate")
variable Mode1D=SrwMode
variable Prec1D=SrwPrec
variable MethTransvEmit=srwUtiGetValN("MethTransvEmit",1,"SrwIntArbMCCreate")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt NumPt,"Number of Macro-Particles"
prompt Mode1D,SrwPMode,popup SrwPOPUPMode
prompt Prec1D,"Relative Prec. for 1D Integration"
prompt MethTransvEmit,"Method for treating Transverse Emittance",popup "by repeated integ. along traj.;by convolution"
Silent 1						|	Computing the Radiation  ...
//PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwIntArbMCCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwIntArbMCCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

srwUtiSetValS("RadName",RadName,"SrwIntArbMCCreate")
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwViewRadCmpnType=RadCmpnType
SrwMagGenTotName=MagName
srwUtiSetValN("NumPt", NumPt, "SrwIntArbMCCreate")
SrwMode=Mode1D
SrwPrec=Prec1D
srwUtiSetValN("MethTransvEmit",MethTransvEmit,"SrwIntArbMCCreate")

variable savePeriodNpt = 100 //to edit

variable elecSigXe2 = $ElecName[20]   // <(x-<x>)^2> [m^2]
variable elecMXXp = $ElecName[21]   // <(x-<x>)(x'-<x'>)>
variable elecSigXpe2 = $ElecName[22]   // <(x'-<x'>)^2>

variable elecSigZe2 = $ElecName[23]   // <(z-<z>)^2> [m^2]
variable elecMZZp = $ElecName[24]   // <(z-<z>)(z'-<z'>)>
variable elecSigZpe2 = $ElecName[25]   // <(z'-<z'>)^2>

variable elecEnGeV = $ElecName[0]
variable elecSigRelE = $ElecName[13]   // relative rms energy spread
variable elecSigE_GeV = elecEnGeV*elecSigRelE

variable r = $ObsName[4]   // observation distance

variable elecSigXeffE2 = elecSigXe2 + r*r*elecSigXpe2 + 2*r*elecMXXp
variable elecSigZeffE2 = elecSigZe2 + r*r*elecSigZpe2 + 2*r*elecMZZp
variable elecSigXeff = sqrt(elecSigXeffE2)
variable elecSigZeff = sqrt(elecSigZeffE2)
variable elecSigXeff_mm = elecSigXeff*1000
variable elecSigZeff_mm = elecSigZeff*1000

variable obsEstart_eV = $ObsName[5]
variable obsEend_eV = $ObsName[6]
variable obsEnp = $ObsName[7]
variable obsXstart_m = $ObsName[8]
variable obsXend_m = $ObsName[9]
variable obsXnp = $ObsName[10]
variable obsZstart_m = $ObsName[11]
variable obsZend_m = $ObsName[12]
variable obsZnp = $ObsName[13]

variable obsEstart_keV = obsEstart_eV*0.001
variable obsEend_keV = obsEend_eV*0.001
variable obsXmid_mm = 0.5*(obsXstart_m + obsXend_m)*1000
variable obsXrange_mm = (obsXend_m - obsXstart_m)*1000
variable obsZmid_mm = 0.5*(obsZstart_m + obsZend_m)*1000
variable obsZrange_mm = (obsZend_m - obsZstart_m)*1000

variable obsXmidAct = obsXmid_mm, obsZmidAct = obsZmid_mm
variable elecEnAct = elecEnGeV

string auxElecName = "ElecAuxMC_ebm"
duplicate/O $ElecName $auxElecName

string auxObsName = "ObsAuxMC"
SrwSmpCreate(auxObsName,r)
auxObsName += "_obs"

string auxWfrName = "WfrAuxMC"
string auxSufIntWork = "Iw"
string auxEndingInt = "_"
variable numDimInt = 0
//to program correctly !
if(obsEnp > 1)
	auxEndingInt += "e"
	numDimInt += 1
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	else
	
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	else
	
	endif
else
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	endif
endif

variable SingE_or_MultiE_intens = 1 //single-e by default
if(MethTransvEmit == 2) //by convolution
	SingE_or_MultiE_intens = 2
endif

string auxIntWorkName = auxWfrName + auxSufIntWork + auxEndingInt
string auxIntResName = RadName + auxEndingInt
string auxFluxDispWorkName = auxWfrName + auxSufIntWork + "_e"
string auxFluxDispResName = RadName + "_e"

variable dispOrNotWork = 2
if(numDimInt == 3)
	dispOrNotWork = 1 //don't display the main wave
	//auxIntDisp += "_e"
endif
//variable eStep, eStart

//if(MethTransvEmit == 2) //use convolution
//	if(IntOrFluxWork == 1)
//		IntOrFluxWork == 2 //multi-e intensity
//	endif
//endif

variable savePtCount = 0
variable iPt = 0
do
	//defining next values for transverse positions of the observation and electron energy
	elecEnAct = elecEnGeV + gnoise(elecSigE_GeV)
	$auxElecName[0] = elecEnAct
	
	if(iPt == 0)
		obsXmidAct = obsXmid_mm
		obsZmidAct = obsZmid_mm
	else
		if(MethTransvEmit == 2) //by convolution
			obsXmidAct = obsXmid_mm
			obsZmidAct = obsZmid_mm
		else
			obsXmidAct = obsXmid_mm - gnoise(elecSigXeff_mm)
			obsZmidAct = obsZmid_mm - gnoise(elecSigZeff_mm)
		endif
	endif	
	
	SrwSmpScanXZE(auxObsName,obsXmidAct,obsXrange_mm,obsXnp,obsZmidAct,obsZrange_mm,obsZnp,obsEstart_keV,obsEend_keV,obsEnp)

	SrwMagPrec(MagName,Mode1D,0.01,Prec1D,10000,1,0,0)
	SrwWfrCreate(auxWfrName,auxElecName,MagName,auxObsName,1,1)
	
	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,RadCmpnType,SingE_or_MultiE_intens,8,1,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	 
	 if(numDimInt == 3)
	 	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,7,5,1,1,obsEstart_keV,obsXmidAct,obsZmidAct,2) //extracting spectral flux for viewing
	 endif
	
	if(iPt == 0)
		SrwUtiGraphWindResize(400,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()
		
		duplicate/O $auxIntWorkName $auxIntResName
		
		if(numDimInt == 1)
			display $auxIntResName
		endif
		if(numDimInt == 2)
			display; AppendImage $auxIntResName
		endif
		if(numDimInt == 3) //extract and show spectrum through finite aperture
			duplicate/O $auxFluxDispWorkName $auxFluxDispResName
			display $auxFluxDispResName
		endif
		
		SrwUtiGraphWindResize(10,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()	
		dispOrNotWork = 1
	else
		$auxIntResName *= iPt
		
		if(numDimInt == 1)
			if(MethTransvEmit == 1)
				$auxIntResName += $auxIntWorkName[p]
			else
				$auxIntResName += $auxIntWorkName(x)
			endif
		else
			if(numDimInt == 2)
				if(MethTransvEmit == 1)
					$auxIntResName += $auxIntWorkName[p][q]
				else
					$auxIntResName += $auxIntWorkName(x)(y)
				endif
			else
				if(MethTransvEmit == 1)
					$auxIntResName += $auxIntWorkName[p][q][r]
				else
					$auxIntResName += $auxIntWorkName(x)(y)(z)
				endif
			endif
		endif
		
		$auxIntResName /= (iPt + 1)
		
		if(numDimInt == 3)
			$auxFluxDispResName *= iPt
			if(MethTransvEmit == 1)
				$auxFluxDispResName += $auxFluxDispWorkName[p]
			else
				$auxFluxDispResName += $auxFluxDispWorkName(x)
			endif
			$auxFluxDispResName /= (iPt + 1)
		endif
	endif
	
	DoUpdate

	savePtCount += 1
	
	if(savePtCount == savePeriodNpt)
		SaveExperiment
		savePtCount = 0
	endif
	
	iPt += 1
while(iPt < NumPt)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculate Intensity of radiation taking into account electron beam emittance and energy spread
//by summing-up contributions from different "macro-electrons"
//Directly integrates over e-beam phase-space volume.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwIntArbMCCreateTest(RadName, ElecName, MagName, ObsName, RadCmpnType, NumPt, Mode1D, Prec1D, MethTransvEmit)
string RadName=srwUtiGetValS("RadName", "IntME", "SrwIntArbMCCreate")
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string ObsName=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
variable RadCmpnType=SrwViewRadCmpnType
variable NumPt=srwUtiGetValN("NumPt", 100, "SrwIntArbMCCreate")
variable Mode1D=SrwMode
variable Prec1D=SrwPrec
variable MethTransvEmit=srwUtiGetValN("MethTransvEmit",1,"SrwIntArbMCCreate")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt NumPt,"Number of Macro-Particles"
prompt Mode1D,SrwPMode,popup SrwPOPUPMode
prompt Prec1D,"Relative Prec. for 1D Integration"
prompt MethTransvEmit,"Method for treating Transverse Emittance",popup "by repeated integ. along traj.;by convolution"
Silent 1						|	Computing the Radiation  ...
//PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwIntArbMCCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwIntArbMCCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

srwUtiSetValS("RadName",RadName,"SrwIntArbMCCreate")
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwViewRadCmpnType=RadCmpnType
SrwMagGenTotName=MagName
srwUtiSetValN("NumPt", NumPt, "SrwIntArbMCCreate")
SrwMode=Mode1D
SrwPrec=Prec1D
srwUtiSetValN("MethTransvEmit",MethTransvEmit,"SrwIntArbMCCreate")

variable savePeriodNpt = 100 //to edit

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

variable elecEnGeV = srwGetElecBeamEnergy(ElecName)
variable elecSigE_GeV = elecEnGeV*elecRelEnSpr

variable elecCur = srwGetElecBeamCurrent(ElecName)
variable elecS0 = srwGetElecBeamLongPos(ElecName)
variable elecX0 = srwGetElecBeamHorPos(ElecName)
variable elecXp0 = srwGetElecBeamHorAng(ElecName)
variable elecZ0 = srwGetElecBeamVertPos(ElecName)
variable elecZp0 = srwGetElecBeamVertAng(ElecName)

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
variable elecX, elecXp, elecZ, elecZp, elecEnAct = elecEnGeV

//variable elecSigXe2 = $ElecName[20]   // <(x-<x>)^2> [m^2]
//variable elecMXXp = $ElecName[21]   // <(x-<x>)(x'-<x'>)>
//variable elecSigXpe2 = $ElecName[22]   // <(x'-<x'>)^2>

//variable elecSigZe2 = $ElecName[23]   // <(z-<z>)^2> [m^2]
//variable elecMZZp = $ElecName[24]   // <(z-<z>)(z'-<z'>)>
//variable elecSigZpe2 = $ElecName[25]   // <(z'-<z'>)^2>

//variable elecEnGeV = $ElecName[0]
//variable elecSigRelE = $ElecName[13]   // relative rms energy spread
//variable elecSigE_GeV = elecEnGeV*elecSigRelE

//variable r = $ObsName[4]   // observation distance
//variable elecSigXeffE2 = elecSigXe2 + r*r*elecSigXpe2 + 2*r*elecMXXp
//variable elecSigZeffE2 = elecSigZe2 + r*r*elecSigZpe2 + 2*r*elecMZZp
//variable elecSigXeff = sqrt(elecSigXeffE2)
//variable elecSigZeff = sqrt(elecSigZeffE2)
//variable elecSigXeff_mm = elecSigXeff*1000
//variable elecSigZeff_mm = elecSigZeff*1000

variable obsEstart_eV = $ObsName[5]
variable obsEend_eV = $ObsName[6]
variable obsEnp = $ObsName[7]
variable obsXstart_m = $ObsName[8]
variable obsXend_m = $ObsName[9]
variable obsXnp = $ObsName[10]
variable obsZstart_m = $ObsName[11]
variable obsZend_m = $ObsName[12]
variable obsZnp = $ObsName[13]

variable obsEstart_keV = obsEstart_eV*0.001
variable obsEend_keV = obsEend_eV*0.001
variable obsXmid_mm = 0.5*(obsXstart_m + obsXend_m)*1000
//variable obsXrange_mm = (obsXend_m - obsXstart_m)*1000
variable obsZmid_mm = 0.5*(obsZstart_m + obsZend_m)*1000
//variable obsZrange_mm = (obsZend_m - obsZstart_m)*1000

//variable obsXmidAct = obsXmid_mm, obsZmidAct = obsZmid_mm

string ElecWorkCore = "ElecAuxMC"
string ElecWorkName = ElecWorkCore + "_ebm"

//string auxElecName = "ElecAuxMC_ebm"
//duplicate/O $ElecName $auxElecName

//string auxObsName = "ObsAuxMC"
//SrwSmpCreate(auxObsName,r)
//auxObsName += "_obs"

string auxWfrName = "WfrAuxMC"
string auxSufIntWork = "Iw"
string auxEndingInt = "_"
variable numDimInt = 0
//to program correctly !
if(obsEnp > 1)
	auxEndingInt += "e"
	numDimInt += 1
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	else
	
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	else
	
	endif
else
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	endif
endif

variable SingE_or_MultiE_intens = 1 //single-e by default
if(MethTransvEmit == 2) //by convolution
	SingE_or_MultiE_intens = 2
endif

string auxIntWorkName = auxWfrName + auxSufIntWork + auxEndingInt
string auxIntResName = RadName + auxEndingInt
string auxFluxDispWorkName = auxWfrName + auxSufIntWork + "_e"
string auxFluxDispResName = RadName + "_e"

variable dispOrNotWork = 2
if(numDimInt == 3)
	dispOrNotWork = 1 //don't display the main wave
	//auxIntDisp += "_e"
endif
//variable eStep, eStart

//if(MethTransvEmit == 2) //use convolution
//	if(IntOrFluxWork == 1)
//		IntOrFluxWork == 2 //multi-e intensity
//	endif
//endif

string nmXcSigLPTau = "wCentersSigmasPropStokesMultiE"
string nmElecInitCond = "wElecInitCondPropStokesMultiE"

make/O $nmXcSigLPTau = {{0,0,0,0,0},{1,1,1,1,1}}
make/O/N=5 $nmElecInitCond
variable initRand = 1, modeRand = 0 // 0- random, 1- LPTau

variable savePtCount = 0
variable iPt = 0
do

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


	//defining next values for transverse positions of the observation and electron energy
	//elecEnAct = elecEnGeV + gnoise(elecSigE_GeV)
	//$auxElecName[0] = elecEnAct
	
	//if(iPt == 0)
	//	obsXmidAct = obsXmid_mm
	//	obsZmidAct = obsZmid_mm
	//else
	//	if(MethTransvEmit == 2) //by convolution
	//		obsXmidAct = obsXmid_mm
	//		obsZmidAct = obsZmid_mm
	//	else
	//		obsXmidAct = obsXmid_mm - gnoise(elecSigXeff_mm)
	//		obsZmidAct = obsZmid_mm - gnoise(elecSigZeff_mm)
	//	endif
	//endif	
	
	elecEnAct = elecEnGeV + elecSigE_GeV*$nmElecInitCond[4]

		//print "i=", iElec, "  Electron Coordinates: x=", elecX, " x'=", elecXp, "z=", elecZ, " z'=", elecZp, "energy=", elecEnAct
	
	SrwElecFilament(ElecWorkCore, elecEnAct, elecCur, elecS0, elecX*1000, elecXp*1000, elecZ*1000, elecZp*1000)

	//SrwSmpScanXZE(auxObsName,obsXmidAct,obsXrange_mm,obsXnp,obsZmidAct,obsZrange_mm,obsZnp,obsEstart_keV,obsEend_keV,obsEnp)

	SrwMagPrec(MagName,Mode1D,0.008,Prec1D,10000,1,0,0)
	//SrwWfrCreate(auxWfrName,auxElecName,MagName,auxObsName,1,1)
	//SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,RadCmpnType,SingE_or_MultiE_intens,8,1,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	
	SrwWfrCreate(auxWfrName,ElecWorkName,MagName,ObsName,1,1)
	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,RadCmpnType,SingE_or_MultiE_intens,8,1,obsEstart_keV,obsXmid_mm,obsZmid_mm,dispOrNotWork)
	
	 if(numDimInt == 3)
	 	//SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,7,5,1,1,obsEstart_keV,obsXmidAct,obsZmidAct,2) //extracting spectral flux for viewing
	 	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,7,5,1,1,obsEstart_keV,obsXmid_mm,obsZmid_mm,2) //extracting spectral flux for viewing
	 endif
	
	if(iPt == 0)
		SrwUtiGraphWindResize(400,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()
		
		duplicate/O $auxIntWorkName $auxIntResName
		
		if(numDimInt == 1)
			display $auxIntResName
		endif
		if(numDimInt == 2)
			display; AppendImage $auxIntResName
		endif
		if(numDimInt == 3) //extract and show spectrum through finite aperture
			duplicate/O $auxFluxDispWorkName $auxFluxDispResName
			display $auxFluxDispResName
		endif
		
		SrwUtiGraphWindResize(10,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()	
		dispOrNotWork = 1
	else
		$auxIntResName *= iPt
		
		if(numDimInt == 1)
			$auxIntResName += $auxIntWorkName(x)
		else
			if(numDimInt == 2)
				$auxIntResName += $auxIntWorkName(x)(y)
			else
				$auxIntResName += $auxIntWorkName(x)(y)(z)
			endif
		endif
		
		$auxIntResName /= (iPt + 1)
		
		if(numDimInt == 3)
			$auxFluxDispResName *= iPt
			$auxFluxDispResName += $auxFluxDispWorkName(x)
			$auxFluxDispResName /= (iPt + 1)
		endif
	endif
	
	DoUpdate

	savePtCount += 1
	
	if(savePtCount == savePeriodNpt)
		SaveExperiment
		savePtCount = 0
	endif
	
	iPt += 1
while(iPt < NumPt)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//TEST: Calculates Flux through finite aperture of radiation from Arbitrary Source taking into account electron beam emittance and energy spread
//by summing-up contributions from different "macro-electrons"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwFluxArbMCCreate(RadName, ElecName, MagName, ObsName, RadCmpnType, NumPt, Mode1D, Prec1D, MethTransvEmit, iCalcMade)
string RadName=srwUtiGetValS("RadName", "IntME", "SrwIntArbMCCreate")
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string ObsName=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
variable RadCmpnType=SrwViewRadCmpnType
variable NumPt=srwUtiGetValN("NumPt", 100, "SrwIntArbMCCreate")
variable Mode1D=SrwMode
variable Prec1D=SrwPrec
variable MethTransvEmit=srwUtiGetValN("MethTransvEmit",1,"SrwIntArbMCCreate")
variable iCalcMade=srwUtiGetValN("iCalcMade",1,"SrwFluxArbMCCreate")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt NumPt,"Number of Macro-Particles"
prompt Mode1D,SrwPMode,popup SrwPOPUPMode
prompt Prec1D,"Relative Prec. for 1D Integration"
prompt MethTransvEmit,"Method for treating Transverse Emittance",popup "by repeated integ. along traj.;by convolution"
prompt iCalcMade,"Number of Part. treated (effective if >=0)"
Silent 1						|	Computing the Radiation  ...
//PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwIntArbMCCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwIntArbMCCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

srwUtiSetValS("RadName",RadName,"SrwIntArbMCCreate")
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwViewRadCmpnType=RadCmpnType
SrwMagGenTotName=MagName
srwUtiSetValN("NumPt", NumPt, "SrwIntArbMCCreate")
SrwMode=Mode1D
SrwPrec=Prec1D
srwUtiSetValN("MethTransvEmit",MethTransvEmit,"SrwIntArbMCCreate")
srwUtiSetValN("iCalcMade",iCalcMade,"SrwIntArbMCCreate")

variable savePeriodNpt = 100 //to edit

variable elecSigXe2 = $ElecName[20]   // <(x-<x>)^2> [m^2]
variable elecMXXp = $ElecName[21]   // <(x-<x>)(x'-<x'>)>
variable elecSigXpe2 = $ElecName[22]   // <(x'-<x'>)^2>

variable elecSigZe2 = $ElecName[23]   // <(z-<z>)^2> [m^2]
variable elecMZZp = $ElecName[24]   // <(z-<z>)(z'-<z'>)>
variable elecSigZpe2 = $ElecName[25]   // <(z'-<z'>)^2>

variable elecEnGeV = $ElecName[0]
variable elecSigRelE = $ElecName[13]   // relative rms energy spread
variable elecSigE_GeV = elecEnGeV*elecSigRelE

variable r = $ObsName[4]   // observation distance

variable elecSigXeffE2 = elecSigXe2 + r*r*elecSigXpe2 + 2*r*elecMXXp
variable elecSigZeffE2 = elecSigZe2 + r*r*elecSigZpe2 + 2*r*elecMZZp
variable elecSigXeff = sqrt(elecSigXeffE2)
variable elecSigZeff = sqrt(elecSigZeffE2)
variable elecSigXeff_mm = elecSigXeff*1000
variable elecSigZeff_mm = elecSigZeff*1000

variable obsEstart_eV = $ObsName[5]
variable obsEend_eV = $ObsName[6]
variable obsEnp = $ObsName[7]
variable obsXstart_m = $ObsName[8]
variable obsXend_m = $ObsName[9]
variable obsXnp = $ObsName[10]
variable obsZstart_m = $ObsName[11]
variable obsZend_m = $ObsName[12]
variable obsZnp = $ObsName[13]

variable obsEstart_keV = obsEstart_eV*0.001
variable obsEend_keV = obsEend_eV*0.001
variable obsXmid_mm = 0.5*(obsXstart_m + obsXend_m)*1000
variable obsXrange_mm = (obsXend_m - obsXstart_m)*1000
variable obsZmid_mm = 0.5*(obsZstart_m + obsZend_m)*1000
variable obsZrange_mm = (obsZend_m - obsZstart_m)*1000

variable obsHalfXrange_mm = 0.5*obsXrange_mm
variable obsHalfZrange_mm = 0.5*obsZrange_mm
variable multInt2Flux = obsXrange_mm*obsZrange_mm

variable obsXmidAct = obsXmid_mm, obsZmidAct = obsZmid_mm
variable elecEnAct = elecEnGeV

string auxElecName = "ElecAuxMC_ebm"
duplicate/O $ElecName $auxElecName

string auxObsName = "ObsAuxMC"
SrwSmpCreate(auxObsName,r)
auxObsName += "_obs"

string auxWfrName = "WfrAuxMC"
string auxSufIntWork = "Iw"
string auxEndingInt = "_"
variable numDimInt = 0
//to program correctly !
if(obsEnp > 1)
	auxEndingInt += "e"
	numDimInt += 1
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	else
	
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	else
	
	endif
else
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	endif
endif

variable SingE_or_MultiE_intens = 1 //single-e by default
if(MethTransvEmit == 2) //by convolution
	SingE_or_MultiE_intens = 2
endif

string auxIntWorkName = auxWfrName + auxSufIntWork + auxEndingInt
string auxIntResName = RadName + auxEndingInt
string auxFluxDispWorkName = auxWfrName + auxSufIntWork + "_e"
string auxFluxDispResName = RadName + "_e"

variable dispOrNotWork = 2
if(numDimInt == 3)
	dispOrNotWork = 1 //don't display the main wave
	//auxIntDisp += "_e"
endif
//variable eStep, eStart

//if(MethTransvEmit == 2) //use convolution
//	if(IntOrFluxWork == 1)
//		IntOrFluxWork == 2 //multi-e intensity
//	endif
//endif

variable savePtCount = 0
variable iPt = 0

if(iCalcMade >= 0) //OC250112
	iPt += iCalcMade
	NumPt += iCalcMade
endif

do
	//defining next values for transverse positions of the observation and electron energy
	elecEnAct = elecEnGeV + gnoise(elecSigE_GeV)
	$auxElecName[0] = elecEnAct
	
	if(MethTransvEmit == 2) //by convolution
		obsXmidAct = obsXmid_mm
		obsZmidAct = obsZmid_mm
	else
		//obsXmidAct = obsXmid_mm - gnoise(elecSigXeff_mm)
		//obsZmidAct = obsZmid_mm - gnoise(elecSigZeff_mm)
		
		obsXmidAct = obsXmid_mm + enoise(obsHalfXrange_mm) - gnoise(elecSigXeff_mm)
		obsZmidAct = obsZmid_mm + enoise(obsHalfZrange_mm) - gnoise(elecSigZeff_mm)
	endif
	SrwSmpScanXZE(auxObsName,obsXmidAct,obsXrange_mm,obsXnp,obsZmidAct,obsZrange_mm,obsZnp,obsEstart_keV,obsEend_keV,obsEnp)

	SrwMagPrec(MagName,Mode1D,0.01,Prec1D,10000,1,0,0)
	SrwWfrCreate(auxWfrName,auxElecName,MagName,auxObsName,1,1)
	
	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,RadCmpnType,SingE_or_MultiE_intens,8,1,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	 
	 if(numDimInt == 3)
	 	SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,7,5,1,1,obsEstart_keV,obsXmidAct,obsZmidAct,2) //extracting spectral flux for viewing
	 endif
	
	if(iPt == 0)
		SrwUtiGraphWindResize(400,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()
		
		duplicate/O $auxIntWorkName $auxIntResName
		
		$auxIntResName *= multInt2Flux
		
		if(numDimInt == 1)
			display $auxIntResName
		endif
		if(numDimInt == 2)
			display; AppendImage $auxIntResName
		endif
		if(numDimInt == 3) //extract and show spectrum through finite aperture
			duplicate/O $auxFluxDispWorkName $auxFluxDispResName
			display $auxFluxDispResName
		endif
		
		SrwUtiGraphWindResize(10,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()	
		dispOrNotWork = 1
	else
		//$auxIntResName *= iPt
		$auxIntResName *= iPt/multInt2Flux
		$auxIntResName += $auxIntWorkName
		//$auxIntResName /= (iPt + 1)
		$auxIntResName *= multInt2Flux/(iPt + 1)
		
		if(numDimInt == 3)
			$auxFluxDispResName *= iPt
			$auxFluxDispResName += $auxFluxDispWorkName
			$auxFluxDispResName /= (iPt + 1)
		endif
	endif
	
	DoUpdate

	savePtCount += 1
	
	if(savePtCount == savePeriodNpt)
		SaveExperiment
		savePtCount = 0
	endif
	
	iPt += 1
while(iPt < NumPt)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//TEST: Calculates Flux through finite aperture of Bending Magnet radiation taking into account electron beam emittance and energy spread
//by summing-up contributions from different "macro-electrons"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//proc SrwFluxBMMCCreate(RadName, ElecName, MagName, ObsName, RadCmpnType, NumPt, Mode1D, Prec1D, MethTransvEmit, iCalcMade)
proc SrwFluxBMMCCreate(RadName, ElecName, MagName, ObsName, RadCmpnType, NumPt, MethTransvEmit, iCalcMade)
string RadName=srwUtiGetValS("RadName", "IntME", "SrwIntArbMCCreate")
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string ObsName=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
variable RadCmpnType=SrwViewRadCmpnType
variable NumPt=srwUtiGetValN("NumPt", 100, "SrwIntArbMCCreate")
//variable Mode1D=SrwMode
//variable Prec1D=SrwPrec
variable MethTransvEmit=srwUtiGetValN("MethTransvEmit",1,"SrwIntArbMCCreate")
variable iCalcMade=srwUtiGetValN("iCalcMade",1,"SrwFluxArbMCCreate")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwMagConstType,";","") //+ Wavelist("*"+SrwMagContainerType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt NumPt,"Number of Macro-Particles"
//prompt Mode1D,SrwPMode,popup SrwPOPUPMode
//prompt Prec1D,"Relative Prec. for 1D Integration"
prompt MethTransvEmit,"Method for treating Transverse Emittance",popup "by repeated integ. along traj.;by convolution"
prompt iCalcMade,"Number of Part. treated (effective if >=0)"
Silent 1						|	Computing the Radiation  ...
//PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwIntArbMCCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwIntArbMCCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

srwUtiSetValS("RadName",RadName,"SrwIntArbMCCreate")
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwViewRadCmpnType=RadCmpnType
SrwMagGenTotName=MagName
srwUtiSetValN("NumPt", NumPt, "SrwIntArbMCCreate")
//SrwMode=Mode1D
//SrwPrec=Prec1D
srwUtiSetValN("MethTransvEmit",MethTransvEmit,"SrwIntArbMCCreate")
srwUtiSetValN("iCalcMade",iCalcMade,"SrwIntArbMCCreate")

variable savePeriodNpt = 100 //to edit

variable elecSigXe2 = $ElecName[20]   // <(x-<x>)^2> [m^2]
variable elecMXXp = $ElecName[21]   // <(x-<x>)(x'-<x'>)>
variable elecSigXpe2 = $ElecName[22]   // <(x'-<x'>)^2>

variable elecSigZe2 = $ElecName[23]   // <(z-<z>)^2> [m^2]
variable elecMZZp = $ElecName[24]   // <(z-<z>)(z'-<z'>)>
variable elecSigZpe2 = $ElecName[25]   // <(z'-<z'>)^2>

variable elecEnGeV = $ElecName[0]
variable elecSigRelE = $ElecName[13]   // relative rms energy spread
variable elecSigE_GeV = elecEnGeV*elecSigRelE

variable r = $ObsName[4]   // observation distance

variable elecSigXeffE2 = elecSigXe2 + r*r*elecSigXpe2 + 2*r*elecMXXp
variable elecSigZeffE2 = elecSigZe2 + r*r*elecSigZpe2 + 2*r*elecMZZp
variable elecSigXeff = sqrt(elecSigXeffE2)
variable elecSigZeff = sqrt(elecSigZeffE2)
variable elecSigXeff_mm = elecSigXeff*1000
variable elecSigZeff_mm = elecSigZeff*1000

variable obsEstart_eV = $ObsName[5]
variable obsEend_eV = $ObsName[6]
variable obsEnp = $ObsName[7]
variable obsXstart_m = $ObsName[8]
variable obsXend_m = $ObsName[9]
variable obsXnp = $ObsName[10]
variable obsZstart_m = $ObsName[11]
variable obsZend_m = $ObsName[12]
variable obsZnp = $ObsName[13]

variable obsEstart_keV = obsEstart_eV*0.001
variable obsEend_keV = obsEend_eV*0.001
variable obsXmid_mm = 0.5*(obsXstart_m + obsXend_m)*1000
variable obsXrange_mm = (obsXend_m - obsXstart_m)*1000
variable obsZmid_mm = 0.5*(obsZstart_m + obsZend_m)*1000
variable obsZrange_mm = (obsZend_m - obsZstart_m)*1000

variable obsHalfXrange_mm = 0.5*obsXrange_mm
variable obsHalfZrange_mm = 0.5*obsZrange_mm
variable multInt2Flux = obsXrange_mm*obsZrange_mm

variable obsXmidAct = obsXmid_mm, obsZmidAct = obsZmid_mm
variable elecEnAct = elecEnGeV

string auxElecName = "ElecAuxMC_ebm"
duplicate/O $ElecName $auxElecName

string auxObsName = "ObsAuxMC"
SrwSmpCreate(auxObsName,r)
auxObsName += "_obs"

string auxWfrName = "WfrAuxMC"
string auxSufIntWork = "Iw"
string auxEndingInt = "_"
variable numDimInt = 0
//to program correctly !
if(obsEnp > 1)
	auxEndingInt += "e"
	numDimInt += 1
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	else
	
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	else
	
	endif
else
	if(obsXnp > 1)
		auxEndingInt += "x"
		numDimInt += 1
	endif
	if(obsZnp > 1)
		auxEndingInt += "z"
		numDimInt += 1
	endif
endif

variable SingE_or_MultiE_intens = 1 //single-e by default
if(MethTransvEmit == 2) //by convolution
	SingE_or_MultiE_intens = 2
endif

string auxIntWorkName = auxWfrName + auxSufIntWork + auxEndingInt
string auxIntResName = RadName + auxEndingInt
string auxFluxDispWorkName = auxWfrName + auxSufIntWork + "_e"
string auxFluxDispResName = RadName + "_e"

variable dispOrNotWork = 2
if(numDimInt == 3)
	dispOrNotWork = 1 //don't display the main wave
	//auxIntDisp += "_e"
endif
//variable eStep, eStart

//if(MethTransvEmit == 2) //use convolution
//	if(IntOrFluxWork == 1)
//		IntOrFluxWork == 2 //multi-e intensity
//	endif
//endif

variable savePtCount = 0
variable iPt = 0

if(iCalcMade >= 0) //OC250112
	iPt += iCalcMade
	NumPt += iCalcMade
endif

do
	//defining next values for transverse positions of the observation and electron energy
	elecEnAct = elecEnGeV + gnoise(elecSigE_GeV)
	$auxElecName[0] = elecEnAct
	
	if(MethTransvEmit == 2) //by convolution
		obsXmidAct = obsXmid_mm
		obsZmidAct = obsZmid_mm
	else
		//obsXmidAct = obsXmid_mm - gnoise(elecSigXeff_mm)
		//obsZmidAct = obsZmid_mm - gnoise(elecSigZeff_mm)
		
		obsXmidAct = obsXmid_mm + enoise(obsHalfXrange_mm) - gnoise(elecSigXeff_mm)
		obsZmidAct = obsZmid_mm + enoise(obsHalfZrange_mm) - gnoise(elecSigZeff_mm)
	endif
	SrwSmpScanXZE(auxObsName,obsXmidAct,obsXrange_mm,obsXnp,obsZmidAct,obsZrange_mm,obsZnp,obsEstart_keV,obsEend_keV,obsEnp)

	//SrwMagPrec(MagName,Mode1D,0.01,Prec1D,10000,1,0,0)
	//SrwWfrCreate(auxWfrName,auxElecName,MagName,auxObsName,1,1)
	
	SrwStoConstCreate(auxWfrName,auxElecName,MagName,auxObsName)
	
	//SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,RadCmpnType,SingE_or_MultiE_intens,8,1,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	SrwSto2IntF(auxWfrName + "_ras",auxSufIntWork,RadCmpnType,1,8,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	 
	 if(numDimInt == 3)
	 	//SrwWfr2Int(auxWfrName + "_rad",auxSufIntWork,7,5,1,1,obsEstart_keV,obsXmidAct,obsZmidAct,2) //extracting spectral flux for viewing
		SrwSto2IntF(auxWfrName + "_ras",auxSufIntWork,7,2,8,obsEstart_keV,obsXmidAct,obsZmidAct,dispOrNotWork)
	 endif
	
	if(iPt == 0)
		SrwUtiGraphWindResize(400,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()
		
		duplicate/O $auxIntWorkName $auxIntResName
		
		$auxIntResName *= multInt2Flux
		
		if(numDimInt == 1)
			display $auxIntResName
		endif
		if(numDimInt == 2)
			display; AppendImage $auxIntResName
		endif
		if(numDimInt == 3) //extract and show spectrum through finite aperture
			duplicate/O $auxFluxDispWorkName $auxFluxDispResName
			display $auxFluxDispResName
		endif
		
		SrwUtiGraphWindResize(10,10,350,200,0,0)
		SrwUtiGraphAddFrameAndGrid()	
		dispOrNotWork = 1
	else
		//$auxIntResName *= iPt
		$auxIntResName *= iPt/multInt2Flux
		$auxIntResName += $auxIntWorkName
		//$auxIntResName /= (iPt + 1)
		$auxIntResName *= multInt2Flux/(iPt + 1)
		
		if(numDimInt == 3)
			$auxFluxDispResName *= iPt
			$auxFluxDispResName += $auxFluxDispWorkName
			$auxFluxDispResName /= (iPt + 1)
		endif
	endif
	
	DoUpdate

	savePtCount += 1
	
	if(savePtCount == savePeriodNpt)
		SaveExperiment
		savePtCount = 0
	endif
	
	iPt += 1
while(iPt < NumPt)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Stokes from Arbitrary Magnetic Field and Thick Electron Beam. 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoArbCreate(RadName, ElecName, MagName, ObsName, IntOrFlux, MethNo, StepOrRelPrec, MaxIter)
string RadName=srwUtiTruncString(SrwElecName+SrwUndName+SrwSmpName, 27)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string ObsName=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
variable IntOrFlux=srwUtiGetValN("IntOrFlux", 1, "SrwStoArbCreate")
variable MethNo=srwUtiGetValN("MethNo", 1, "SrwStoArbCreate")
variable StepOrRelPrec=srwUtiGetValN("StepOrRelPrec", 0.01, "SrwStoArbCreate")
variable MaxIter=srwUtiGetValN("MaxIter", 1000, "SrwStoArbCreate")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt IntOrFlux,"Characteristics to Calculate",popup "Spectral Flux Through a Slit;Spectral Flux per Unit Surface"
prompt MethNo,"Method of Integration",popup "Manual;Automatic;Test"
prompt StepOrRelPrec,"Integ. Step [m] (man.) or Rel. Prec. (auto)"
prompt MaxIter,"Maximal Number of Iterations (auto)"
Silent 1						|	Computing the Radiation  ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwStoArbCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwStoArbCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

// Validation of parameters
//...

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]

SrwMagGenTotName=MagName

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

srwUtiSetValN("IntOrFlux", IntOrFlux, "SrwStoArbCreate")
srwUtiSetValN("MethNo", MethNo, "SrwStoArbCreate")
srwUtiSetValN("StepOrRelPrec", StepOrRelPrec, "SrwStoArbCreate")
srwUtiSetValN("MaxIter", MaxIter, "SrwStoArbCreate")

// Preparing data for C function
Make/D/O/N=6 waveprec
waveprec[0]=IntOrFlux // What to calculate
waveprec[1]=MethNo // Method No
waveprec[2]=StepOrRelPrec // Step or precision
waveprec[3]=MaxIter // Max. number of iterations

SrwStoPrep(ElecName,MagName,Obsname,RadName,IntOrFlux)
RadName += SrwStoType
SrwRadGenTotName=RadName

srStokesArb($ElecName, $MagName, $ObsName, waveprec, $RadName)
KillWaves/Z  waveprec
end

//+++++++++++++++++++++++++++++++++++++++
//
//Auxiliary function: Computes flux of single-electron emission for given electron energy
//to be eventually used to accelerate SrwWfrEmitPropStokesMultiE calculation 
//
//+++++++++++++++++++++++++++++++++++++++
function srwStoArbFluxVsElecEn(enElect, nmElec, nmMag, nmObs, intMeth, relPrec, polCmpn)
variable enElect, intMeth, relPrec, polCmpn
string nmElec, nmMag, nmObs

variable dummyIntStep = 0.01
if(intMeth == 1)
	dummyIntStep = relPrec
endif
DoUpdate

srwSetElecBeamEnergy(nmElec, enElect)
string str2exe = "SrwMagPrec(\""+nmMag+"\","+num2str(intMeth)+","+num2str(dummyIntStep)+","+num2str(relPrec)+",10000,1,0,0);"
str2exe += "SrwWfrCreate(\"auxWfrStoArbFluxVsElecEn\",\""+nmElec+"\",\""+nmMag+"\",\""+nmObs+"\",1,1);"
str2exe += "SrwWfr2Int(\"auxWfrStoArbFluxVsElecEn_rad\",\"F\","+num2str(polCmpn)+",5,1,1,1,0,0,1)"
	//print str2exe
execute/Q str2exe
wave auxWfrStoArbFluxVsElecEnF_e
return auxWfrStoArbFluxVsElecEnF_e[0]
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Dedicated to compute flux through finite aperture at a peak of harmonic at variable K from arbitrary magnetic field and thick electron beam.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwStoArbFluxAtVarK(FluxName,ElecName,MagName,ObsName,Polar,Harm,Prec,RelKRange,NumK,Disp)
string FluxName=srwUtiGetValS("FluxName","FluxVsK","SrwStoArbFluxAtVarK")
string ElecName=srwUtiGetValS("SrwElecName","Elec","")+SrwElecType
string MagName=srwUtiGetValS("SrwMagName","Mag","")+SrwFieldType
string ObsName=srwUtiGetValS("SrwSmpName","Obs","")+SrwSmpType
variable Harm=srwUtiGetValN("Harm",1,"SrwStoArbFluxAtVarK")
variable Polar=srwUtiGetValN("Polar",1,"SrwStoArbFluxAtVarK")
variable Prec=srwUtiGetValN("Prec",0.01,"SrwStoArbFluxAtVarK")
variable RelKRange=srwUtiGetValN("RelKRange",0.1,"SrwStoArbFluxAtVarK")
variable NumK=srwUtiGetValN("NumK",30,"SrwStoArbFluxAtVarK")
variable Disp=srwUtiGetValN("Disp",1,"SrwStoArbFluxAtVarK")
prompt FluxName,"Name for the Flux vs K structure"
prompt ElecName,"Electron Beam structure",popup wavelist("*"+SrwElecType,";","")
prompt MagName,"Magnetic Field at Max. K structure",popup wavelist("*"+SrwFieldType,";","")
prompt ObsName,"Radiation Sampling structure",popup wavelist("*"+SrwSmpType,";","")
prompt Harm,"Spectrum Harmonic number"
prompt Polar,"Polarization Component",popup SrwPOPUPPolar+";Total"
prompt Prec,"Relative Precision"
prompt RelKRange,"Rel. Range of K values to find Max."
prompt NumK,"Number of K values to find Max."
prompt Disp,"New Display ?",popup SrwPOPUPViewDisplay
Silent 1						|	Computing spectral flux at harmonic vs K ...
PauseUpdate

if((strlen(FluxName)==0) %| (cmpstr(FluxName,"_none_")==0))
	abort "The value of the parameter \"Name for the Flux vs K structure\" is not defined"
endif
if(strlen(FluxName) > 28)
	abort "The string parameter \"Name for the Flux vs K structure\" is too long"
endif
if((strlen(ElecName)==0) %| (cmpstr(ElecName,"_none_")==0))
	abort "The value of the parameter \"Electron Beam structure\" is not defined"
endif
if((strlen(MagName)==0) %| (cmpstr(MagName,"_none_")==0))
	abort "The value of the parameter \"Magnetic Field at Max. K structure\" is not defined"
endif
if((strlen(ObsName)==0) %| (cmpstr(ObsName,"_none_")==0))
	abort "The value of the parameter \"Radiation Sampling structure\" is not defined"
endif

srwUtiSetValS("FluxName",FluxName,"SrwStoArbFluxAtVarK")
srwUtiSetValS("SrwElecName",ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1],"")
srwUtiSetValS("SrwMagName",MagName[0,strlen(MagName)-strlen(SrwFieldType)-1],"")
srwUtiSetValS("SrwSmpName",ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1],"")
srwUtiSetValN("Harm",Harm,"SrwStoArbFluxAtVarK")
srwUtiSetValN("Polar",Polar,"SrwStoArbFluxAtVarK")
srwUtiSetValN("Prec",Prec,"SrwStoArbFluxAtVarK")
srwUtiSetValN("RelKRange",RelKRange,"SrwStoArbFluxAtVarK")
srwUtiSetValN("NumK",NumK,"SrwStoArbFluxAtVarK")
srwUtiSetValN("Disp",Disp,"SrwStoArbFluxAtVarK")

variable xc = 0.5*(srwGetSmpHorPosStart(ObsName) + srwGetSmpHorPosEnd(ObsName))
variable xr =  srwGetSmpHorPosEnd(ObsName) - srwGetSmpHorPosStart(ObsName)
variable nx = srwGetSmpHorPosNp(ObsName)
variable zc = 0.5*(srwGetSmpVertPosStart(ObsName) + srwGetSmpVertPosEnd(ObsName))
variable zr =  srwGetSmpVertPosEnd(ObsName) - srwGetSmpVertPosStart(ObsName)
variable nz = srwGetSmpVertPosNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
variable ne = srwGetSmpPhotEnNp(ObsName)
variable eStep = (eEnd - eStart)/(ne - 1)

string NameFieldModulCoefs = FluxName + "MC"
make/O/N=(ne) $FluxName, $NameFieldModulCoefs
SetScale/I x eStart,eEnd,"eV", $FluxName, $NameFieldModulCoefs
$FluxName = 0
$NameFieldModulCoefs = 0
variable/G OptFieldModulCoef

SrwSmpCreate("ObsAux",srwGetSmpLongPos(ObsName))

if(Disp == 2) 
	display $NameFieldModulCoefs; SrwUtiGraphAddFrameAndGrid()
	display $FluxName; SrwUtiGraphAddFrameAndGrid()
	SrwUtiGraphWindResize(150,10,350,200,0,0)
	DoUpdate
endif

variable ie = 0
variable CurPhEn = eStart
variable MaxIterElecEnSpread = 1 //3

do
	SrwSmpScanXZE("ObsAux_obs",xc*1000,xr*1000,nx,zc*1000,zr*1000,nz,CurPhEn*0.001,CurPhEn*0.001,1)
	$FluxName[ie] = SrwStoArbFluxAtFixEnVarK(ElecName,MagName,"ObsAux_obs",Polar,Harm,Prec,MaxIterElecEnSpread,RelKRange,NumK,1)
	$NameFieldModulCoefs[ie] = OptFieldModulCoef
	DoUpdate
	
	CurPhEn += eStep
	ie += 1
while(ie < ne)

end 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Dedicated to compute flux through finite aperture at a peak of harmonic at variable K from arbitrary magnetic field and thick electron beam.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SrwStoArbFluxAtFixEnVarK(ElecName,MagName,ObsName,Polar,Harm,Prec,MaxIter,RelRangeK,NumK,AbsMaxFieldModCoef)
string ElecName, MagName, ObsName
variable Harm, Polar, Prec, MaxIter, RelRangeK, NumK, AbsMaxFieldModCoef

string/G SrwFieldType

//Convert arbitrary field to periodic to estimate fundamental photon energy
SVAR SrwFieldType
string PerFldName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1] + "P"

string StrToExe = "SrwMagArb2Per(\"" + PerFldName + "\",\"" + MagName + "\",0.05,10,1000)"
execute StrToExe

variable Kx1 = 0, Kz1 = 0
string HarmName1 = PerFldName + "1" + "_fha"
string HarmName2 = PerFldName + "2" + "_fha"
string HarmName3 = PerFldName + "3" + "_fha"
string HarmName4 = PerFldName + "4" + "_fha"
if(exists(HarmName1) == 1)
	wave wHarmName1 = $HarmName1
	if(wHarmName1[0] == 1)
		if(wHarmName1[1] == 1)
			Kz1 = wHarmName1[2]
		else 
			Kx1 = wHarmName1[2]			
		endif
	endif
endif
if(exists(HarmName2) == 1)
	wave wHarmName2 = $HarmName2
	if(wHarmName2[0] == 1)
		if(wHarmName2[1] == 1)
			Kz1 = wHarmName2[2]
		else 
			Kx1 = wHarmName2[2]			
		endif
	endif
endif
if(exists(HarmName3) == 1)
	wave wHarmName3 = $HarmName3
	if(wHarmName3[0] == 1)
		if(wHarmName3[1] == 1)
			Kz1 = wHarmName3[2]
		else 
			Kx1 = wHarmName3[2]			
		endif
	endif
endif
if(exists(HarmName4) == 1)
	wave wHarmName4 = $HarmName4
	if(wHarmName4[0] == 1)
		if(wHarmName4[1] == 1)
			Kz1 = wHarmName4[2]
		else 
			Kx1 = wHarmName4[2]			
		endif
	endif
endif

SVAR SrwUndType

variable EffOrigK = sqrt(Kx1*Kx1 + Kz1*Kz1)
string TotPerFldName = PerFldName + SrwUndType
wave/T wTotPerFld = $TotPerFldName
variable UndPer = str2num(wTotPerFld[0])
variable UndNper = round(str2num(wTotPerFld[1])/UndPer)

wave wElec = $ElecName
variable ElecEn = wElec[0]

variable a1 = 9.5
variable e1max = a1*ElecEn*ElecEn/UndPer
variable enmax = e1max*Harm

variable PhotEn = srwGetSmpPhotEnStart(ObsName)
if(PhotEn > enmax)
	return 0
endif

//if(exists("PrevK") == 0)
variable/G PrevK = 0
//else
//	variable/G PrevK
//endif

variable EstimK = 0
variable EstimValE1 = PhotEn/Harm + PhotEn/UndNper
//if(PrevK == 0)
EstimK = sqrt(2*(a1*ElecEn*ElecEn/(EstimValE1*UndPer) - 1))
//endif

variable DeltaE1 = RelRangeK*EstimValE1
variable DeltaK = DeltaE1*abs(-a1*ElecEn*ElecEn/(EstimValE1*EstimValE1*UndPer)/sqrt(2*(a1*ElecEn*ElecEn/(EstimValE1*UndPer) - 1)))
//RelRangeK = abs(-a1*ElecEn*ElecEn/(EstimValE1*EstimValE1*UndPer)/sqrt(2*(a1*ElecEn*ElecEn/(EstimValE1*UndPer) - 1)))
//variable DeltaK = RelRangeK*EstimK

variable Kmin = EstimK - DeltaK, Kmax = EstimK + DeltaK
if(Kmin < 0.0001)
	Kmin = 0.0001
endif

variable FieldModulCoefMin = Kmin/EffOrigK, FieldModulCoefMax = Kmax/EffOrigK
if(FieldModulCoefMax > AbsMaxFieldModCoef)
	FieldModulCoefMax = AbsMaxFieldModCoef
endif
if(FieldModulCoefMin > AbsMaxFieldModCoef)
	FieldModulCoefMin = AbsMaxFieldModCoef
endif

variable FieldModulCoefStep = (FieldModulCoefMax - FieldModulCoefMin)/(NumK - 1)
variable CurFieldModulCoef = FieldModulCoefMin
if(FieldModulCoefStep <= 0)
	NumK = 1
endif
variable ik = 0

string OrigMagFldNameCore = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
string NewMagFldNameCore = OrigMagFldNameCore + "D"
string NewMagFldNameTot = NewMagFldNameCore + "_mag"
StrToExe = "SrwMagDupl(\"" + MagName + "\",\"" + NewMagFldNameCore + "\")"
execute StrToExe

string OrigNameBx = OrigMagFldNameCore + "BX_fld", OrigNameBz = OrigMagFldNameCore + "BZ_fld"
string NameBx = NewMagFldNameCore + "BX_fld", NameBz = NewMagFldNameCore + "BZ_fld"
wave wOrigBx = $OrigNameBx, wOrigBz = $OrigNameBz, wBx = $NameBx, wBz = $NameBz

variable MaxFlux = 0
variable/G OptFieldModulCoef = -1
string AuxFluxName = "StokesAuxI_e"

do
	wBx = CurFieldModulCoef*wOrigBx[p]
	wBz = CurFieldModulCoef*wOrigBz[p]
	
	StrToExe = "DoUpdate; SrwStoArbCreate(\"StokesAux\",\"" + ElecName + "\",\"" + NewMagFldNameTot + "\",\"" + ObsName + "\",2,3," + num2str(Prec) + "," + num2str(MaxIter) + ")"
	execute StrToExe
	StrToExe = "SrwSto2IntF(\"StokesAux_ras\",\"I\"," + num2str(Polar) + ",2,1," + num2str(PhotEn) + ",0,0,1)"
	execute StrToExe
	
	wave wInt = $AuxFluxName
	if(MaxFlux < wInt[0])
		MaxFlux = wInt[0]
		OptFieldModulCoef = CurFieldModulCoef
	endif

	CurFieldModulCoef += FieldModulCoefStep
	ik += 1
while(ik < NumK)

return MaxFlux
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Dedicated to compute flux through finite aperture at a peak of harmonic at variable K from arbitrary magnetic field and thick electron beam.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SrwStoArbFluxAtHarmPeak(ElecName,MagName,ObsName,Polar,Harm,Prec,RelPhotEnRange,NumPhEn)
string ElecName, MagName, ObsName
variable Harm, Polar, Prec, RelPhotEnRange, NumPhEn

SVAR SrwFieldType
//Convert arbitrary field to periodic to estimate fundamental photon energy
string PerFldName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1] + "P"

string StrToExe = "SrwMagArb2Per(\"" + PerFldName + "\",\"" + MagName + "\",0.05,10,1000)"
execute StrToExe

variable Kx1 = 0, Kz1 = 0
string HarmName1 = PerFldName + "1" + "_fha"
string HarmName2 = PerFldName + "2" + "_fha"
string HarmName3 = PerFldName + "3" + "_fha"
string HarmName4 = PerFldName + "4" + "_fha"
if(exists(HarmName1) == 1)
	wave wHarmName1 = $HarmName1
	if(wHarmName1[0] == 1)
		if(wHarmName1[1] == 1)
			Kz1 = wHarmName1[2]
		else 
			Kx1 = wHarmName1[2]			
		endif
	endif
endif
if(exists(HarmName2) == 1)
	wave wHarmName2 = $HarmName2
	if(wHarmName2[0] == 1)
		if(wHarmName2[1] == 1)
			Kz1 = wHarmName2[2]
		else 
			Kx1 = wHarmName2[2]			
		endif
	endif
endif
if(exists(HarmName3) == 1)
	wave wHarmName3 = $HarmName3
	if(wHarmName3[0] == 1)
		if(wHarmName3[1] == 1)
			Kz1 = wHarmName3[2]
		else 
			Kx1 = wHarmName3[2]			
		endif
	endif
endif
if(exists(HarmName4) == 1)
	wave wHarmName4 = $HarmName4
	if(wHarmName4[0] == 1)
		if(wHarmName4[1] == 1)
			Kz1 = wHarmName4[2]
		else 
			Kx1 = wHarmName4[2]			
		endif
	endif
endif

SVAR SrwUndType

string TotPerFldName = PerFldName + SrwUndType
wave/T wTotPerFld = $TotPerFldName
variable UndPer = str2num(wTotPerFld[0])
wave wElec = $ElecName
variable ElecEn = wElec[0]
variable e1 = 9.5*ElecEn*ElecEn/(UndPer*(1 + 0.5*(Kx1*Kx1 + Kz1*Kz1)))

variable eDelta = RelPhotEnRange*e1
variable ec = e1*Harm
variable eMin = ec - eDelta, eMax = ec // + eDelta
variable eStep = (eMax - eMin)/(NumPhEn - 1)
variable ie = 0

StrToExe = "SrwSmpCreate(\"ObsAux\"," + num2str(srwGetSmpLongPos(ObsName)) + ")"
execute StrToExe

variable xc = 0.5*(srwGetSmpHorPosStart(ObsName) + srwGetSmpHorPosEnd(ObsName))
variable xr =  srwGetSmpHorPosEnd(ObsName) - srwGetSmpHorPosStart(ObsName)
variable nx = srwGetSmpHorPosNp(ObsName)
variable zc = 0.5*(srwGetSmpVertPosStart(ObsName) + srwGetSmpVertPosEnd(ObsName))
variable zr =  srwGetSmpVertPosEnd(ObsName) - srwGetSmpVertPosStart(ObsName)
variable nz = srwGetSmpVertPosNp(ObsName)
variable PhotEn = eMin
variable MaxFlux = 0
string AuxFluxName = "StokesAuxI_e"
do
	StrToExe = "SrwSmpScanXZE(\"ObsAux_obs\"," + num2str(xc*1000) + "," + num2str(xr*1000) + "," + num2str(nx) + "," + num2str(zc*1000) + "," + num2str(zr*1000) + "," + num2str(nz) + "," + num2str(PhotEn*0.001) + "," + num2str(PhotEn*0.001) + ",1)"
	execute StrToExe
	StrToExe = "SrwStoArbCreate(\"StokesAux\",\"" + ElecName + "\",\"" + MagName + "\",\"ObsAux_obs\",2,3," + num2str(Prec) + ")"
	execute StrToExe
	StrToExe = "SrwSto2IntF(\"StokesAux_ras\",\"I\"," + num2str(Polar) + ",2,1," + num2str(PhotEn) + ",0,0,1)"
	execute StrToExe
	
	wave wInt = $AuxFluxName
	if(MaxFlux < wInt[0])
		MaxFlux = wInt[0]
	endif

	PhotEn += eStep
	ie += 1
while(ie < NumPhEn)
return MaxFlux
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Processes single-electron intensity calculated vs photon energy, horizontal and vertical position in order
//to take into account e-beam energy spread
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwIntTreatElecEnSprAndEmit(nmIntensMD, nmElec, nmMag, obsDist_m, RadCmpnType, Mode1D, Prec1D, toTreat)
string nmIntensMD=srwUtiGetValS("nmIntensMD", "IntMD", "SrwIntTreatElecEnSprAndEmit")
string nmElec=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string nmMag=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
variable obsDist_m=srwUtiGetValN("obsDist_m", 20, "SrwIntTreatElecEnSprAndEmit")
variable RadCmpnType=SrwViewRadCmpnType
variable Mode1D=SrwMode
variable Prec1D=SrwPrec
variable toTreat=srwUtiGetValN("toTreat", 1, "SrwIntTreatElecEnSprAndEmit")
prompt nmIntensMD,"Single-Electron Intensity (vs e, x, z) to treat",popup Wavelist("*",";","TEXT:0,DIMS:3")
prompt nmElec,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt nmMag,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt obsDist_m,"Observation Distance [m] for which the Intensity was calculated"
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt Mode1D,SrwPMode,popup SrwPOPUPMode
prompt Prec1D,"Relative Prec. for 1D Integration"
prompt toTreat,"To take into account:",popup "Energy Spread only;Transverse Emittance only;Energy Spread and Transverse Emittance"
Silent 1						|	Computing Radiation  ...
PauseUpdate

srwUtiSetValS("nmIntensMD",nmIntensMD,"SrwIntTreatElecEnSprAndEmit")
SrwElecName=nmElec[0,strlen(nmElec)-strlen(SrwElecType)-1]
SrwUndName=nmMag[0,strlen(nmMag)-strlen(SrwUndType)-1]
srwUtiSetValN("obsDist_m", obsDist_m, "SrwIntTreatElecEnSprAndEmit")
SrwViewRadCmpnType=RadCmpnType
SrwMagGenTotName=nmMag
SrwMode=Mode1D
SrwPrec=Prec1D
srwUtiSetValN("toTreat", toTreat, "SrwIntTreatElecEnSprAndEmit")

variable elecEn0_GeV = srwGetElecBeamEnergy(nmElec)
variable elecSigRelE = srwGetElecBeamRelEnSprRMS(nmElec)  // relative rms energy spread
variable elecSigE_GeV = elecEn0_GeV*elecSigRelE
variable elecEnOff_GeV = elecEn0_GeV + elecSigE_GeV
variable elecX0_m = srwGetElecBeamHorPos(nmElec)
variable elecXp0_r = srwGetElecBeamHorAng(nmElec)
variable elecZ0_m = srwGetElecBeamVertPos(nmElec)
variable elecZp0_r = srwGetElecBeamVertAng(nmElec)
variable elecS0_m = srwGetElecBeamLongPos(nmElec)
variable elecCurrent_A = srwGetElecBeamCurrent(nmElec)

variable elecSigX = srwGetElecBeamHorSizeRMS(nmElec) //$ElecName[20]   // <(x-<x>)^2> [m^2]
variable elecMXXp = srwGetElecBeamHorMixedMom(nmElec) //$ElecName[21]   // <(x-<x>)(x'-<x'>)>
variable elecSigXp = srwGetElecBeamHorDivergRMS(nmElec) //$ElecName[22]   // <(x'-<x'>)^2>
variable elecSigZ = srwGetElecBeamVertSizeRMS(nmElec) //$ElecName[23]   // <(z-<z>)^2> [m^2]
variable elecMZZp = srwGetElecBeamVertMixedMom(nmElec) //$ElecName[24]   // <(z-<z>)(z'-<z'>)>
variable elecSigZp = srwGetElecBeamVertDivergRMS(nmElec) //$ElecName[25]   // <(z'-<z'>)^2>

variable elecXproj = elecX0_m + obsDist_m*elecXp0_r
variable elecZproj = elecZ0_m + obsDist_m*elecZp0_r
variable elecSigXe2 = elecSigX*elecSigX
variable elecSigXpe2 = elecSigXp*elecSigXp
variable elecSigZe2 = elecSigZ*elecSigZ
variable elecSigZpe2 = elecSigZp*elecSigZp

variable obsEstart_eV = dimoffset($nmIntensMD, 0)
variable obsEstep_eV = dimdelta($nmIntensMD, 0)
variable obsEnp = dimsize($nmIntensMD, 0)
variable obsEend_eV = obsEstart_eV + (obsEnp - 1)*obsEstep_eV
variable obsEstart_keV = obsEstart_eV*0.001, obsEend_keV = obsEend_eV*0.001

variable obsXstart_m = dimoffset($nmIntensMD, 1)
variable obsXstep_m = dimdelta($nmIntensMD, 1)
variable obsXnp = dimsize($nmIntensMD, 1)
variable obsXrange_m = (obsXnp - 1)*obsXstep_m
variable obsXcen_m = obsXstart_m + 0.5*obsXrange_m
variable obsX0_m = obsXcen_m - elecXproj
variable obsXoff_m = obsX0_m + 0.5*obsXrange_m
variable ix0 = round((obsX0_m - obsXstart_m)/obsXstep_m)
if(ix0 >= obsXnp)
	ix0 = obsXnp - 1
endif
obsX0_m = obsXstart_m + ix0*obsXstep_m //ensure that the test point traps exactly onto a mesh point
variable ixOff = round((obsXoff_m - obsXstart_m)/obsXstep_m)
if(ixOff >= obsXnp)
	ixOff = obsXnp - 1
endif
obsXoff_m = obsXstart_m + ixOff*obsXstep_m //ensure that the test point traps exactly onto a mesh point
variable obsX0_mm = obsX0_m*1000, obsXoff_mm = obsXoff_m*1000

variable obsZstart_m = dimoffset($nmIntensMD, 2)
variable obsZstep_m = dimdelta($nmIntensMD, 2)
variable obsZnp = dimsize($nmIntensMD, 2)
variable obsZrange_m = (obsZnp - 1)*obsZstep_m
variable obsZcen_m = obsZstart_m + 0.5*obsZrange_m
variable obsZ0_m = obsZcen_m - elecZproj
variable obsZoff_m = obsZ0_m + 0.5*obsZrange_m
variable iz0 = round((obsZ0_m - obsZstart_m)/obsZstep_m)
if(iz0 >= obsZnp)
	iz0 = obsZnp - 1
endif
obsZ0_m = obsZstart_m + iz0*obsZstep_m //ensure that the test point traps exactly onto a mesh point
variable izOff = round((obsZoff_m - obsZstart_m)/obsZstep_m)
if(izOff >= obsZnp)
	izOff = obsZnp - 1
endif
obsZoff_m = obsZstart_m + izOff*obsZstep_m //ensure that the test point traps exactly onto a mesh point
variable obsZ0_mm = obsZ0_m*1000, obsZoff_mm = obsZoff_m*1000

variable obsOffAngX = (obsXoff_m - elecXproj)/obsDist_m
variable obsOffAngZ = (obsZoff_m - elecZproj)/obsDist_m
variable sqSumOffAng = obsOffAngX*obsOffAngX + obsOffAngZ*obsOffAngZ

//setup e-beam with different electron energy
SrwElecFilament(SrwElecName, elecEnOff_GeV, elecCurrent_A, elecS0_m, elecX0_m*1000, elecXp0_r*1000, elecZ0_m*1000, elecZp0_r*1000)

string nmAuxObs = "AuxObsTreatElecEnSpr", nmAuxWfr = "AuxWfrTreatElecEnSpr"
string sufAuxInt0 = "I0", sufAuxIntOff = "IOff"
string nmInt0 = nmAuxWfr + sufAuxInt0 + "_e", nmIntOff = nmAuxWfr + sufAuxIntOff + "_e"
string nmInt0Orig = "or" + nmInt0, nmIntOffOrig = "or" + nmIntOff

SrwSmpCreate(nmAuxObs, obsDist_m)

//setup parameters and calculate on-axis spectrum with different electron energy
SrwSmpScanXZE(nmAuxObs + "_obs", obsX0_mm, 1, 1, obsZ0_mm, 1, 1, obsEstart_keV, obsEend_keV, obsEnp)
SrwMagPrec(nmMag, Mode1D, Prec1D, Prec1D, 10000, 1, 0, 0)
SrwWfrCreate(nmAuxWfr, nmElec, nmMag, nmAuxObs + "_obs", 1, 1)
SrwWfr2Int(nmAuxWfr + "_rad", sufAuxInt0, RadCmpnType, 1, 1, 1, obsEstart_keV, obsX0_mm, obsZ0_mm, 1) //extract intensity without displaying
duplicate/O $nmInt0 $nmInt0Orig
$nmInt0Orig = $nmIntensMD[p][ix0][iz0]

//setup parameters and calculate off-axis spectrum with different electron energy
SrwSmpScanXZE(nmAuxObs + "_obs", obsXoff_mm, 1, 1, obsZoff_mm, 1, 1, obsEstart_keV, obsEend_keV, obsEnp)
SrwMagPrec(nmMag, Mode1D, Prec1D, Prec1D, 10000, 1, 0, 0)
SrwWfrCreate(nmAuxWfr, nmElec, nmMag, nmAuxObs + "_obs", 1, 1)
SrwWfr2Int(nmAuxWfr + "_rad", sufAuxIntOff, RadCmpnType, 1, 1, 1, obsEstart_keV, obsXoff_mm, obsZoff_mm, 1) //extract intensity without displaying
duplicate/O $nmIntOff $nmIntOffOrig
$nmIntOffOrig = $nmIntensMD[p][ixOff][izOff]

variable elecRelEnDif = (elecEnOff_GeV - elecEn0_GeV)/elecEn0_GeV
variable halfRangeArgExtCoef = 4*abs(elecRelEnDif) //to tune
variable absTolArgExtCoef = halfRangeArgExtCoef*1e-04
variable a0 = srwUtiFindArgExtCoef($nmInt0Orig, $nmInt0, halfRangeArgExtCoef, absTolArgExtCoef)
variable aOff = srwUtiFindArgExtCoef($nmIntOffOrig, $nmIntOff, halfRangeArgExtCoef, absTolArgExtCoef)
variable b = sqSumOffAng*aOff/(a0 - aOff)
variable p0 = -a0/elecRelEnDif
variable SigE_p0 = p0*elecSigRelE
variable SigE_p

variable numSigPrec = 6 //to tune
variable iz = 0, ix
obsZoff_m = obsZstart_m

if((toTreat == 1) %| (toTreat == 3))
	do
		ix = 0
		obsXoff_m = obsXstart_m
		do
			$nmIntOffOrig = $nmIntensMD[p][ix][iz]
		
			obsOffAngX = (obsXoff_m - elecXproj)/obsDist_m
			obsOffAngZ = (obsZoff_m - elecZproj)/obsDist_m
			sqSumOffAng = obsOffAngX*obsOffAngX + obsOffAngZ*obsOffAngZ
			SigE_p = SigE_p0*(1 - sqSumOffAng/(b + sqSumOffAng))
			SrwUtiConvWaveWithGausLinVar1D($nmIntOffOrig, SigE_p, numSigPrec, $nmIntOff)
			$nmIntensMD[][ix][iz] = $nmIntOff[p]
	
			obsXoff_m += obsXstep_m
			ix += 1
		while(ix < obsXnp)
		obsZoff_m += obsZstep_m
		iz += 1
	while(iz < obsZnp)
endif
//setting back e-beam with original electron energy
SrwElecFilament(SrwElecName, elecEn0_GeV, elecCurrent_A, elecS0_m, elecX0_m*1000, elecXp0_r*1000, elecZ0_m*1000, elecZp0_r*1000)

//treatting transverse emittance by convolution
variable elecSigXeffE2 = elecSigXe2 + obsDist_m*obsDist_m*elecSigXpe2 + 2*obsDist_m*elecMXXp
variable elecSigZeffE2 = elecSigZe2 + obsDist_m*obsDist_m*elecSigZpe2 + 2*obsDist_m*elecMZZp
if((elecSigXeffE2 <= 0) %& (elecSigZeffE2 <= 0))
	return
endif
if((toTreat != 2) %& (toTreat != 3))
	return
endif

variable elecSigXeff = sqrt(elecSigXeffE2)
variable elecSigZeff = sqrt(elecSigZeffE2)

string nmTreatElecEnSprAndEmit = "wTreatElecEnSprAndEmit_xz"
make/O/N=(obsXnp, obsZnp) $nmTreatElecEnSprAndEmit
SetScale/P x obsXstart_m, obsXstep_m, "m", $nmTreatElecEnSprAndEmit
SetScale/P y obsZstart_m, obsZstep_m, "m", $nmTreatElecEnSprAndEmit

variable ie = 0
do
	$nmTreatElecEnSprAndEmit = $nmIntensMD[ie][p][q]
		//DoUpdate
	SrwUtiConvWaveWithGaus2D(nmTreatElecEnSprAndEmit, elecSigXeff, elecSigZeff)
		//DoUpdate
	$nmIntensMD[ie] = $nmTreatElecEnSprAndEmit[q][r]
	ie += 1
while(ie < obsEnp)

KillWaves/Z $nmTreatElecEnSprAndEmit
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculates "Gap Spectrum", i.e. loops on:
//Simulation of Magnetic Field vs Gap
//Calculation of Intensity on 3D mesh;
//Extraction of Multi-Electron Intensity (without Energy Spread);
//Extraction of Spectral Flux (integrated within fixed aperture);
//Integration of Spectral Flux with a Transmission/Reflection Spectrum
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUndParamScanSpec(nmResSpec, nmElec, nmMag0, undParam0, nmBPeakVsUndParam, nmUndParams, nmObs, nmSpecTr, RadCmpnType, Prec1D)
string nmResSpec = srwUtiGetValS("nmResSpec", "nmResSpec", "SrwUndParamScanSpec")
string nmElec = srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
string nmMag0 = srwUtiGetValS("SrwMagGenTotName", "Mag", "")
variable undParam0 = srwUtiGetValN("undParam0", 1, "SrwUndParamScanSpec")
string nmBPeakVsUndParam = srwUtiGetValS("nmBPeakVsUndParam", "nmBPeakVsUndParam", "SrwUndParamScanSpec")
string nmUndParams = srwUtiGetValS("nmUndParams", "nmUndParams", "SrwUndParamScanSpec")
string nmObs = srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
string nmSpecTr = srwUtiGetValS("nmSpecTr", "nmSpecTr", "SrwUndParamScanSpec")
variable RadCmpnType = SrwViewRadCmpnType
variable Prec1D=SrwPrec
prompt nmResSpec, "Scaled 1D wave vs Und. Param, for the Spectrum"
prompt nmElec, SrwPElecName1, popup Wavelist("*"+SrwElecType,";","")
prompt nmMag0, "Basic Magnetic Field structure", popup Wavelist("*"+SrwFieldType,";","")
prompt undParam0, "Undulator Parameter for the Basic Field"
prompt nmBPeakVsUndParam, "Dependence of Field vs Und. Param.", popup Wavelist("*",";","TEXT:1,DIMS:1")
prompt nmUndParams, "Sequence of Und. Param. values", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt nmObs, SrwPSmpName2, popup Wavelist("*"+SrwSmpType,";","")
prompt nmSpecTr, "Transmission / Reflection Spectrum", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt Prec1D, "Relative Prec. for 1D Integration"
Silent 1						|	Computing Radiation  ...
//PauseUpdate

srwUtiSetValS("nmResSpec", nmResSpec, "SrwUndParamScanSpec")
SrwElecName=nmElec[0,strlen(nmElec)-strlen(SrwElecType)-1]
SrwMagGenTotName=nmMag0
SrwUndName=nmMag0[0,strlen(nmMag0)-strlen(SrwUndType)-1]
srwUtiSetValN("undParam0", undParam0, "SrwUndParamScanSpec")
srwUtiSetValS("nmBPeakVsUndParam", nmBPeakVsUndParam, "SrwUndParamScanSpec")
srwUtiSetValS("nmUndParams", nmUndParams, "SrwUndParamScanSpec")
SrwSmpName=nmObs[0,strlen(nmObs)-strlen(SrwSmpType)-1]
srwUtiSetValS("nmSpecTr", nmSpecTr, "SrwUndParamScanSpec")
SrwViewRadCmpnType=RadCmpnType
SrwPrec=Prec1D

variable npUndParam = dimsize($nmUndParams, 0)
make/O/N=(npUndParam) $nmResSpec

display $nmResSpec vs $nmUndParams
SrwUtiGraphAddFrameAndGrid()

string nmMagAux = "AuxMagUndParamScanSpec"
string nmObsAux = "AuxObsUndParamScanSpec"
string nmWfrAux = "AuxWfrUndParamScanSpec"
string sufIntAux = "Im"
string nmIntAux = nmWfrAux + sufIntAux + "_exz"
string nmSpecAux = "AuxSpecUndParamScan_e"
string nmFluxAux = "AuxFluxUndParamScan"

variable xMin = 1000*srwGetSmpHorPosStart(nmObs)
variable xMax = 1000*srwGetSmpHorPosEnd(nmObs)
variable zMin = 1000*srwGetSmpVertPosStart(nmObs)
variable zMax = 1000*srwGetSmpVertPosEnd(nmObs)

variable numExtraSigProj = 3 //to tune
SrwSmpExtForElecBeam(nmObsAux, nmObs, nmElec, numExtraSigProj)
//to set photon energy range at least to that of reflection spectrum?

variable iParam = 0, undParam
do
	undParam = $nmUndParams[iParam]
	SrwMagDuplScale(nmMagAux, undParam, nmMag0, undParam0, $nmBPeakVsUndParam[0], $nmBPeakVsUndParam[1])

	SrwMagPrec(nmMagAux + SrwFieldType, 3, Prec1D, Prec1D, 10000, 1, 0, 0)
	SrwWfrCreate(nmWfrAux, nmElec, nmMagAux + SrwFieldType, nmObsAux + SrwSmpType, 1, 1)
	
	SrwWfr2Int(nmWfrAux + SrwRadType, sufIntAux, RadCmpnType, 2, 7, 1, 1, 0, 0, 1) //extract multi-e intensity on 3D mesh, without treating energy spread
	SrwRadIntensIntegXZ(nmSpecAux, nmIntAux, xMin, xMax, zMin, zMax, 1)
	SrwRadIntensIntegVsPhotEnSpec(nmFluxAux, nmSpecAux, nmSpecTr)
	$nmResSpec[iParam] = $nmFluxAux[0]
		DoUpdate
	
	iParam += 1
while(iParam < npUndParam)
end