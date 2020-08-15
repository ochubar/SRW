
//+++++++++++++++++++++++++++++++++++++++
//
//Visualize Panel
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwVisualizeDialog()

SrwVisualizeCreateBufVars()
DoWindow/K SrwVisualizePanel
SrwVisualizePanel()

end

//+++++++++++++++++++++++++++++++++++++++
Window SrwVisualizePanel() : Panel

PauseUpdate; Silent 1		// building window...
//NewPanel /W=(410,65,758,434) as SrwPVisualizeTitle
NewPanel /W=(410,65,818,434) as SrwPVisualizeTitle
SetDrawLayer UserBack

SrwVisualizeDrawAllContr()

EndMacro

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeDrawAllContr()

variable VertOffset=0
variable RectHorSize=398

String AllRadWaves=Wavelist("*"+SrwRadType,";","")+Wavelist("*"+SrwStoType,";","")+Wavelist("*"+SrwPowType,";","")
Variable ItemNo=sRWaveListItemNo(AllRadWaves, ";", SrwRadGenTotName)
if(ItemNo == 0)
	ItemNo=1
	SrwRadGenTotName=sRWaveListItemName(AllRadWaves, ";", 1)
endif

Variable RadWaveExists=1
if(cmpstr(SrwRadGenTotName,"") == 0)
	RadWaveExists=0
endif

// Computed Radiation structure
//DrawRect 10,8+VertOffset,338,75+VertOffset
DrawRect 10,8+VertOffset,RectHorSize,75+VertOffset

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,25+VertOffset,SrwPRadAllKind
PopupMenu popup0Visualize,pos={25,31+VertOffset},size={255,21}
PopupMenu popup0Visualize,value=#"Wavelist(\"*\"+SrwRadType,\";\",\"\")+Wavelist(\"*\"+SrwStoType,\";\",\"\")+Wavelist(\"*\"+SrwPowType,\";\",\"\")",mode=ItemNo,proc=SrwVisualizeSelectPopupProc

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,70+VertOffset,"Type:"

variable/G SrwVisualizeFluxFromIntStokes = 0
variable StokesFluxCanBeExtr = 0

SrwVisualizeShowEnOnlyBuf = 0
SrwVisualizeRadTypeBuf=1

variable npxaux=0, npzaux=0

if(RadWaveExists != 0)
	string ChosenRadType=SrwRadGenTotName[strlen(SrwRadGenTotName)-strlen(SrwRadType),strlen(SrwRadGenTotName)-1]
	string ChosenStructTitle=SrwRadTitleStokes
	if(cmpstr(ChosenRadType,SrwRadType)==0)
		ChosenStructTitle=SrwRadTitleElecField
		SrwVisualizeRadTypeBuf=1 // Electric Field
		
		SrwVisualizeFluxCanBeExtrBuf = 1
		string AuxExName = $SrwRadGenTotName[0]
		if((dimsize($AuxExName,1) <= 1) %| (dimsize($AuxExName,2) <= 1))
			SrwVisualizeFluxCanBeExtrBuf = 0
		endif
	endif
	if(cmpstr(ChosenRadType,SrwStoType)==0)
		ChosenStructTitle=SrwRadTitleStokes
		SrwVisualizeRadTypeBuf=2 // Stokes
		
		string DataUn=WaveUnits($SrwRadGenTotName,-1)
		if(cmpstr(DataUn, SrwPUnitSpAngFlux) == 0)
			StokesFluxCanBeExtr = 1
		endif
		
		//if(cmpstr(DataUn, SrwPUnitSpAngFluxPerUnSurf) == 0)
		npxaux=dimsize($SrwRadGenTotName,2)
		npzaux=dimsize($SrwRadGenTotName,3)
		if((npxaux > 1) %& (npzaux > 1))
			StokesFluxCanBeExtr = 1
		endif
		//endif
	endif
	if(cmpstr(ChosenRadType,SrwPowType)==0)
		ChosenStructTitle=SrwRadTitlePowDens
		SrwVisualizeRadTypeBuf=3 // Power Density
	endif
	SetDrawEnv fname= "default",fsize= 12,fstyle=0
	DrawText 80,70+VertOffset,ChosenStructTitle
endif

// Component
//DrawRect 10,80+VertOffset,338,165+VertOffset
DrawRect 10,80+VertOffset,RectHorSize,165+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,97+VertOffset,"Component"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,118+VertOffset,"Polarization"
PopupMenu popup1Visualize,pos={180,98+VertOffset},size={255,21}
if((SrwVisualizeRadTypeBuf==1) %| (SrwVisualizeRadTypeBuf==2))
	PopupMenu popup1Visualize,value=#"SrwPOPUPPolar+\";Total\"",mode=SrwViewRadCmpnType,proc=SrwVisualizePolarPopupProc
endif
if(SrwVisualizeRadTypeBuf==3)
	PopupMenu popup1Visualize,value="Total",mode=1
	SrwViewRadCmpnType=7
endif

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,139+VertOffset,"Single-e or Multi-e"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,160+VertOffset,"Extract ..."

variable PrevFluxCanBeExtrBuf = srwUtiGetValN("SrwVisualizeFluxCanBeExtr_G", 0, "")

if(SrwVisualizeRadTypeBuf==1) // El. Field

	PopupMenu popup2Visualize,pos={180,119+VertOffset},size={255,21}
	PopupMenu popup2Visualize,value="Single-e;Multi-e",mode=SrwVisualizeSingleOrMultiEBuf,proc=SrwVisualizeSorMPopupProc

	PopupMenu popup5Visualize,pos={180,140+VertOffset},size={255,21}
	
	if(SrwViewRadCmpnType != 7)
		if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e
		
			if(SrwVisualizeFluxCanBeExtrBuf == 0)
				if((PrevFluxCanBeExtrBuf > 0) %& (SrwVisualizeWhatToExtrBuf >= 2))
					//SrwVisualizeWhatToExtrBuf -= 1
					SrwVisualizeWhatToExtrBuf = 1
				endif
				
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 4, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Polarization Rate I1/(I1+I2);Phase;Re(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 5, 1)
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 6, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Pol. Rate I1/(I1+I2);Pol. Rate (I1-I2)/(I1+I2);Phase;Re(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Pol. Rate I1/(I1+I2);Pol. Rate (I1-I2)/(I1+I2);Phase;Re(E);Im(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			else
				if((PrevFluxCanBeExtrBuf <= 0) %& (SrwVisualizeWhatToExtrBuf >= 2))
					SrwVisualizeWhatToExtrBuf += 1
				endif
				
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 5, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Spectral Flux;Polariz. Rate I1/(I1+I2);Phase;Re(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 8, 1)
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 9, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux;Pol. Rate I1/(I1+I2) from Int.;Pol. Rate (I1-I2)/(I1+I2) from Int.;Pol. Rate F1/(F1+F2) from Flux;Pol. Rate (F1-F2)/(F1+F2) from Flux;Phase;Re(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux;Pol. Rate I1/(I1+I2) from Int.;Pol. Rate (I1-I2)/(I1+I2) from Int.;Pol. Rate F1/(F1+F2) from Flux;Pol. Rate (F1-F2)/(F1+F2) from Flux;Phase;Re(E);Im(E)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			endif
			
		else // Multi-e
			//if(SrwVisualizeWhatToExtrBuf>2)
			if(SrwVisualizeWhatToExtrBuf > 6)
				SrwVisualizeWhatToExtrBuf=1
			endif
			
			if(SrwVisualizeFluxCanBeExtrBuf == 0)
				if((PrevFluxCanBeExtrBuf > 0) %& (SrwVisualizeWhatToExtrBuf >= 2))
					//SrwVisualizeWhatToExtrBuf -= 1
					SrwVisualizeWhatToExtrBuf = 1
				endif
				
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 2, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Polarization Rate I1/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 3, 1)
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Pol. Rate I1/(I1+I2);Pol. Rate (I1-I2)/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			else
				if((PrevFluxCanBeExtrBuf <= 0) %& (SrwVisualizeWhatToExtrBuf >= 2))
					SrwVisualizeWhatToExtrBuf += 1
				endif
				
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 3, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Spectral Flux;Polarization Rate I1/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 6, 1)
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux;Pol. Rate I1/(I1+I2) from Int.;Pol. Rate (I1-I2)/(I1+I2) from Int.;Pol. Rate F1/(F1+F2) from Flux;Pol. Rate (F1-F2)/(F1+F2) from Flux",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			endif
		endif
	else
		if(SrwVisualizeFluxCanBeExtrBuf == 0)
			SrwVisualizeWhatToExtrBuf=1
			
			SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 1, 1)
			PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity)",mode=SrwVisualizeWhatToExtrRes
		else
			if((PrevFluxCanBeExtrBuf <= 0) %& (SrwVisualizeWhatToExtrBuf >= 2))
				SrwVisualizeWhatToExtrBuf = 1
			endif
			
			SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 2, 1)
			PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
		endif
	endif
	
	if((SrwVisualizeFluxCanBeExtrBuf > 0) %& (SrwVisualizeWhatToExtrBuf == 2)) //flux chosen
		SrwVisualizeShowEnOnlyBuf = 1
	else
		if(SrwVisualizeFluxCanBeExtrBuf > 0)
			if((SrwVisualizeWhatToExtrBuf == 5) %| (SrwVisualizeWhatToExtrBuf == 6))
				SrwVisualizeShowEnOnlyBuf = 1
			endif
		else
			SrwVisualizeShowEnOnlyBuf = 0
		endif
	endif
	
endif

variable/G SrwVisualizeStokesIntAndFlux=0
if(SrwVisualizeRadTypeBuf==2) // Stokes

	PopupMenu popup2Visualize,pos={180,119+VertOffset},size={255,21}
	PopupMenu popup2Visualize,value="Multi-e",mode=1
	SrwVisualizeSingleOrMultiEBuf=2

	PopupMenu popup5Visualize,pos={180,140+VertOffset},size={255,21}

	if(StokesFluxCanBeExtr == 0)
		if(SrwViewRadCmpnType != 7)
		
			//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 2, 1)
			//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Polarization Rate I1/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 3, 1)
			PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Pol. Rate I1/(I1+I2);Pol. Rate (I1-I2)/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
		else
			SrwVisualizeWhatToExtrBuf=1
			SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 1, 1)
			PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity)",mode=SrwVisualizeWhatToExtrRes
		endif
	else
		if((npxaux <= 1) %| (npzaux <= 1))
			if(SrwViewRadCmpnType != 7)
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 2, 1)
				//PopupMenu popup5Visualize,value="Spectral Flux;Polarization Rate I1/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 3, 1)
				PopupMenu popup5Visualize,value="Spectral Flux;Pol. Rate F1/(F1+F2);Pol. Rate (F1-F2)/(F1+F2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			else
				SrwVisualizeWhatToExtrBuf=1
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 1, 1)
				PopupMenu popup5Visualize,value="Spectral Flux",mode=SrwVisualizeWhatToExtrRes
			endif
		else //stokes on multi-point transverse grid
			if(SrwViewRadCmpnType != 7)
				//SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 3, 1)
				//PopupMenu popup5Visualize,value="Spect. Flux / Surface;Spectral Flux;Polariz. Rate I1/(I1+I2)",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 6, 1)
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux;Pol. Rate I1/(I1+I2) from Int.;Pol. Rate (I1-I2)/(I1+I2) from Int.;Pol. Rate F1/(F1+F2) from Flux;Pol. Rate (F1-F2)/(F1+F2) from Flux",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			else
				SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 2, 1)
				PopupMenu popup5Visualize,value="Spect. Flux / Surface (Intensity);Spectral Flux",mode=SrwVisualizeWhatToExtrRes,proc=SrwVisualizeWhatExtrPopupProc
			endif
			SrwVisualizeStokesIntAndFlux = 1
			
			if(SrwVisualizeWhatToExtrRes == 2) //flux chosen
				SrwVisualizeShowEnOnlyBuf = 1
				SrwVisualizeFluxFromIntStokes = 1
			else
				if((SrwVisualizeWhatToExtrRes == 5) %| (SrwVisualizeWhatToExtrRes == 6)) //polariz. rate via flux chosen
					SrwVisualizeShowEnOnlyBuf = 1
					SrwVisualizeFluxFromIntStokes = 1
				else
					SrwVisualizeShowEnOnlyBuf = 0
				endif
			endif
		endif
	endif
	
endif
if(SrwVisualizeRadTypeBuf==3) // Power

	PopupMenu popup2Visualize,pos={180,119+VertOffset},size={255,21}
	PopupMenu popup2Visualize,value="Multi-e",mode=1
	SrwVisualizeSingleOrMultiEBuf=2

	PopupMenu popup5Visualize,pos={180,140+VertOffset},size={255,21}
	
	SrwVisualizeWhatToExtrBuf=1
	SrwVisualizeWhatToExtrRes = srwUtiNumOrDefault(SrwVisualizeWhatToExtrBuf, 1, 1, 1)
	PopupMenu popup5Visualize,value="Power / Surface",mode=SrwVisualizeWhatToExtrRes
endif

// As a function of
//DrawRect 10,170+VertOffset,338,265+VertOffset
DrawRect 10,170+VertOffset,RectHorSize,265+VertOffset

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,187+VertOffset,SrwPViewPlotType
PopupMenu popup3Visualize,pos={180,175+VertOffset},size={255,21}

if(SrwVisualizeShowEnOnlyBuf == 1)
	PopupMenu popup3Visualize,value="Energy",mode=1,proc=SrwVisualizeAsFuncPopupProc
	SrwVisualizeAsFuncPopupLogics(1)
else
	if((SrwVisualizeRadTypeBuf==1) %| (SrwVisualizeRadTypeBuf==2))
		PopupMenu popup3Visualize,value=#"SrwPOPUPViewPlotType",mode=SrwViewPlotType,proc=SrwVisualizeAsFuncPopupProc
		SrwVisualizeAsFuncPopupLogics(SrwViewPlotType)
	endif
	if(SrwVisualizeRadTypeBuf==3)
		PopupMenu popup3Visualize,value=#"SrwPOPUPViewPlotTypeXZ",mode=SrwViewPlotTypeXZ,proc=SrwVisualizeAsFuncPopupProc
		SrwVisualizeAsFuncPopupLogics(SrwViewPlotTypeXZ)
	endif
endif

Variable OneValue
String NumRad
if(SrwVisualizeShowEnBuf != 0)

	OneValue=0
	if(RadWaveExists != 0)
		if(SrwVisualizeRadTypeBuf == 1) // Electric Field
			NumRad=$SrwRadGenTotName[0]
			if(DimSize($NumRad, 0) == 1)
				SrwViewE = 0.001*DimOffset($NumRad, 0)
				OneValue=1
			endif
		endif
		if(SrwVisualizeRadTypeBuf == 2) // Stokes
			if(DimSize($SrwRadGenTotName, 1) == 1)
				SrwViewE = 0.001*DimOffset($SrwRadGenTotName, 1)
				OneValue=1
			endif
		endif
	endif

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	//DrawText 30,218+VertOffset,SrwPViewE
	DrawText 30,218+VertOffset,"Const. Photon Energy [keV]"
	DrawText 30,229+VertOffset,"           (or Time [fs])"
	
	if(OneValue==0)
		SetVariable setvar1Visualize,pos={180,200+VertOffset},size={120,14},title=" ",fSize=14
		SetVariable setvar1Visualize,limits={-Inf,Inf,1},value=SrwViewE
	else
		SetDrawEnv fname= "default",fsize=14,fstyle=0
		DrawText 183,218+VertOffset, num2str(SrwViewE)
	endif
endif
if(SrwVisualizeShowXBuf != 0)

	OneValue=0
	if(RadWaveExists != 0)
		if(SrwVisualizeRadTypeBuf == 1) // Electric Field
			NumRad=$SrwRadGenTotName[0]
			if(DimSize($NumRad, 1) == 1)
				SrwViewX = 1000.*(DimOffset($NumRad, 1) + 0.5*DimDelta($NumRad, 1))
				OneValue=1
			endif
			if(SrwVisualizeFluxCanBeExtrBuf > 0)
				if(SrwVisualizeWhatToExtrBuf == 2) //flux chosen
					OneValue=1
				endif
				if((SrwVisualizeWhatToExtrBuf == 5) %| (SrwVisualizeWhatToExtrBuf == 6)) //pol. rate via flux chosen
					OneValue=1
				endif
			endif
		endif
		if(SrwVisualizeRadTypeBuf == 2) // Stokes
			if(DimSize($SrwRadGenTotName, 2) == 1)
				SrwViewX = 1000.*(DimOffset($SrwRadGenTotName, 2) + 0.5*DimDelta($SrwRadGenTotName, 2))
				OneValue=1
				//print SrwViewX
			endif
			if(SrwVisualizeShowEnOnlyBuf == 1)
				OneValue=1
			endif
		endif
		if(SrwVisualizeRadTypeBuf == 3) // Power Density
			if(DimSize($SrwRadGenTotName, 0) == 1)
				SrwViewX = 1000.*(DimOffset($SrwRadGenTotName, 0) + 0.5*DimDelta($SrwRadGenTotName, 0))
				OneValue=1
			endif
		endif
	endif

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 30,240+VertOffset,SrwPViewX
	if(OneValue==0)
		SetVariable setvar2Visualize,pos={180,220+VertOffset},size={120,14},title=" ",fSize=14
		SetVariable setvar2Visualize,limits={-Inf,Inf,1},value=SrwViewX
	else
		SetDrawEnv fname= "default",fsize=14,fstyle=0
		DrawText 183,238+VertOffset, num2str(SrwViewX)
	endif
endif
if(SrwVisualizeShowZBuf != 0)

	OneValue=0
	if(RadWaveExists != 0)
		if(SrwVisualizeRadTypeBuf == 1) // Electric Field
			NumRad=$SrwRadGenTotName[0]
			if(DimSize($NumRad, 2) == 1)
				SrwViewZ = 1000.*(DimOffset($NumRad, 2) + 0.5*DimDelta($NumRad, 2))
				OneValue=1
			endif
			if(SrwVisualizeFluxCanBeExtrBuf > 0)
				if(SrwVisualizeWhatToExtrBuf == 2) //flux chosen
					OneValue=1
				endif
				if((SrwVisualizeWhatToExtrBuf == 5) %| (SrwVisualizeWhatToExtrBuf == 6)) //pol. rate via flux chosen
					OneValue=1
				endif				
			endif
		endif
		if(SrwVisualizeRadTypeBuf == 2) // Stokes
			if(DimSize($SrwRadGenTotName, 3) == 1)
				SrwViewZ = 1000.*(DimOffset($SrwRadGenTotName, 3) + 0.5*DimDelta($SrwRadGenTotName, 3))
				OneValue=1
			endif
			if(SrwVisualizeShowEnOnlyBuf == 1)
				OneValue=1
			endif
		endif
		if(SrwVisualizeRadTypeBuf == 3) // Power Density
			if(DimSize($SrwRadGenTotName, 1) == 1)
				SrwViewZ = 1000.*(DimOffset($SrwRadGenTotName, 1) + 0.5*DimDelta($SrwRadGenTotName, 1))
				OneValue=1
			endif
		endif
	endif

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 30,260+VertOffset,SrwPViewZ
	if(OneValue==0)
		SetVariable setvar3Visualize,pos={180,240+VertOffset},size={120,14},title=" ",fSize=14
		SetVariable setvar3Visualize,limits={-Inf,Inf,1},value=SrwViewZ
	else
		SetDrawEnv fname= "default",fsize=14,fstyle=0
		DrawText 183,258+VertOffset, num2str(SrwViewZ)
	endif
endif

// Suffix & Display
//DrawRect 10,270+VertOffset,338,322+VertOffset
DrawRect 10,270+VertOffset,RectHorSize,322+VertOffset

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,287+VertOffset,SrwPViewSuffixExtract
SetVariable setvar0Visualize,pos={25,295+VertOffset},size={130,17},title=" ",fSize=14
SetVariable setvar0Visualize,limits={-Inf,Inf,1},value=SrwSuffixExtract

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 175,287+VertOffset,SrwPViewDisplay
PopupMenu popup4Visualize,pos={180,295+VertOffset},size={255,21}
PopupMenu popup4Visualize,value=#"SrwPOPUPViewDisplay",mode=SrwViewDisplay,proc=SrwVisualizeDispPopupProc

// OK-Cancel-Help
Button button0Visualize,pos={55,337+VertOffset},size={70,20},proc=SrwVisualizeCancelButtonProc,title=SrwPPanelCancelButton
Button button1Visualize,pos={170,337+VertOffset},size={70,20},proc=SrwVisualizeOKButtonProc,title=SrwPPanelOKButton
Button button2Visualize,pos={285,337+VertOffset},size={70,20},proc=SrwVisualizeHelpButtonProc,title=SrwPPanelHelpButton

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeKillAllContr()

KillControl popup0Visualize
KillControl popup1Visualize
KillControl popup2Visualize
KillControl popup3Visualize
KillControl popup4Visualize
KillControl popup5Visualize
KillControl setvar0Visualize
KillControl setvar1Visualize
KillControl setvar2Visualize
KillControl setvar3Visualize

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeSelectPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwRadGenTotName = popStr

srwUtiSetValN("SrwVisualizeFluxCanBeExtr_G", SrwVisualizeFluxCanBeExtrBuf, "")

//SrwVisualizeShowEnOnlyBuf = 0
//if((SrwVisualizeFluxCanBeExtrBuf > 0) %& (SrwVisualizeWhatToExtrBuf == 2)) //flux chosen
//	SrwVisualizeWhatToExtrBuf = 1
//endif

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwVisualizeKillAllContr()
SrwVisualizeDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeAsFuncPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwVisualizeAsFuncPopupLogics(popNum)

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwVisualizeKillAllContr()
SrwVisualizeDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeAsFuncPopupLogics(popNum)
Variable popNum		// which item is currently selected (1-based)

Variable InPopNum=popNum
if((SrwVisualizeRadTypeBuf==1) %| (SrwVisualizeRadTypeBuf==2))
// "Energy;Hor.;Ver.;Hor.+Ver.;En.+Hor.;En.+Ver.;En.+Hor.+Ver.;Auto"
	SrwViewPlotType=InPopNum
	if(popNum==1)
		SrwVisualizeShowEnBuf=0; SrwVisualizeShowXBuf=1; SrwVisualizeShowZBuf=1
	endif
	if(popNum==2)
		SrwVisualizeShowEnBuf=1; SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=1
	endif
	if(popNum==3)
		SrwVisualizeShowEnBuf=1; SrwVisualizeShowXBuf=1; SrwVisualizeShowZBuf=0
	endif
	if(popNum==4)
		SrwVisualizeShowEnBuf=1; SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=0
	endif
	if(popNum==5)
		SrwVisualizeShowEnBuf=0; SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=1
	endif
	if(popNum==6)
		SrwVisualizeShowEnBuf=0; SrwVisualizeShowXBuf=1; SrwVisualizeShowZBuf=0
	endif
	if(popNum==7)
		SrwVisualizeShowEnBuf=0; SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=0
	endif
	if(popNum==8)
		SrwVisualizeShowEnBuf=0; SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=0
	endif
endif
if(SrwVisualizeRadTypeBuf==3)
// "Hor.;Ver.;Hor.+Ver.;Auto"
	SrwViewPlotTypeXZ=InPopNum
	SrwVisualizeShowEnBuf=0
	if(popNum==1)
		SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=1
	endif
	if(popNum==2)
		SrwVisualizeShowXBuf=1; SrwVisualizeShowZBuf=0
	endif
	if(popNum==3)
		SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=0
	endif
	if(popNum==4)
		SrwVisualizeShowXBuf=0; SrwVisualizeShowZBuf=0
	endif
endif
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeDispPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwViewDisplay=popNum

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizePolarPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwViewRadCmpnType=popNum

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwVisualizeKillAllContr()
SrwVisualizeDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeSorMPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwVisualizeSingleOrMultiEBuf=popNum

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwVisualizeKillAllContr()
SrwVisualizeDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeWhatExtrPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwVisualizeWhatToExtrBuf=popNum

srwUtiSetValN("SrwVisualizeFluxCanBeExtr_G", SrwVisualizeFluxCanBeExtrBuf, "")

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwVisualizeKillAllContr()
SrwVisualizeDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeCreateBufVars()

Variable/G SrwVisualizeRadTypeBuf
Variable/G SrwVisualizeShowEnBuf=1, SrwVisualizeShowXBuf=1, SrwVisualizeShowZBuf=1
Variable/G SrwVisualizeShowEnOnlyBuf=0

Variable/G SrwVisualizeFluxCanBeExtrBuf = 1 //srwUtiGetValN("SrwVisualizeFluxCanBeExtr_G", 1, "")
Variable/G SrwVisualizeSingleOrMultiEBuf = srwUtiGetValN("SrwVisualizeSingleOrMultiE_G", 1, "")
Variable/G SrwVisualizeWhatToExtrBuf = srwUtiGetValN("SrwVisualizeWhatToExtr_G", 1, "")
Variable/G SrwVisualizeWhatToExtrRes = SrwVisualizeWhatToExtrBuf

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeKillBufVars()

KillVariables/Z SrwVisualizeRadTypeBuf
KillVariables/Z  SrwVisualizeShowEnBuf, SrwVisualizeShowXBuf, SrwVisualizeShowZBuf
KillVariables/Z SrwVisualizeSingleOrMultiEBuf, SrwVisualizeWhatToExtrBuf

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeCancelButtonProc(ctrlName) : ButtonControl
String ctrlName;

SrwVisualizeDestroy()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeDestroy()

SrwVisualizeKillAllContr()
DoWindow/K SrwVisualizePanel
SrwVisualizeKillBufVars()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeValidateGlobVars()

if(cmpstr(SrwRadGenTotName,"")==0)
	SrwVisualizeDestroy()
	Abort SrwPAlertRadiationNeeded
endif

// Do not allow mutual intensity if only one transverse position
if((SrwVisualizeRadTypeBuf==1) %& (SrwVisualizeSingleOrMultiEBuf==2)) // Electric Field and Multi-e
	String NumWaveName=$SrwRadGenTotName[0]
	Variable nx=DimSize($NumWaveName,1)
	Variable nz=DimSize($NumWaveName,2)
	if((nx<=2) %| (nz<=2))
		Abort SrwPAlertVisMultiE
	endif
endif

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeOKButtonProc(ctrlName) : ButtonControl
String ctrlName

SrwVisualizeWhatToExtrBuf = SrwVisualizeWhatToExtrRes

SrwVisualizeValidateGlobVars()

variable CallWfr2Int = 0
variable CallWfr2PolRate = 0
variable CallSto2Int = 0
variable CallSto2PolRate = 0
variable CallPow2Int = 0

variable PolRateType = 0

variable/G SrwVisualizeStokesIntAndFlux

// Resolve logic Sigle-e/Multi-e, what to extract (to modify later)
if(SrwVisualizeRadTypeBuf==1) // Electric Field
	if(SrwVisualizeFluxCanBeExtrBuf == 0)
		//if(SrwVisualizeWhatToExtrBuf==2)
		if((SrwVisualizeWhatToExtrBuf==2) %| (SrwVisualizeWhatToExtrBuf==3))
			CallWfr2PolRate = 1
			PolRateType = SrwVisualizeWhatToExtrBuf - 1
			
			if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e Intensity
				SrwCmpnNo=1
			else // Multi-e Intensity
				SrwCmpnNo=2
			endif
		else	
			//if((SrwVisualizeWhatToExtrBuf==3) %| (SrwVisualizeWhatToExtrBuf==4)) // Phase or Re(E)
			if((SrwVisualizeWhatToExtrBuf==4) %| (SrwVisualizeWhatToExtrBuf==5)) // Phase or Re(E)
				//SrwCmpnNo=SrwVisualizeWhatToExtrBuf
				SrwCmpnNo=SrwVisualizeWhatToExtrBuf - 1
			else
				if(SrwVisualizeWhatToExtrBuf==6) //Im(E)
					SrwCmpnNo=7
				else
					if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e Intensity
						SrwCmpnNo=1
					else // Multi-e Intensity
						SrwCmpnNo=2
					endif
				endif
			endif
			CallWfr2Int = 1
		endif
	else
		//if(SrwVisualizeWhatToExtrBuf==3)
		if((SrwVisualizeWhatToExtrBuf==3) %| (SrwVisualizeWhatToExtrBuf==4)) //Pol. Rate via Int.
			CallWfr2PolRate = 1
			PolRateType = SrwVisualizeWhatToExtrBuf - 2
			
			if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e Intensity
				SrwCmpnNo=1
			else // Multi-e Intensity
				SrwCmpnNo=2
			endif

		else
			if((SrwVisualizeWhatToExtrBuf==5) %| (SrwVisualizeWhatToExtrBuf==6)) //Pol. Rate via Flux F1/(F1+F2);Pol. Rate via Flux (F1-F2)/(F1+F2)
				CallWfr2PolRate = 1
				
				if(SrwVisualizeWhatToExtrBuf==5)
					PolRateType = 1
				else
					PolRateType = 2
				endif
				
				if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e
					SrwCmpnNo=5 // Single-e Flux
				else // Multi-e
					SrwCmpnNo=6 // Multi-e Flux
				endif
			else
				//if((SrwVisualizeWhatToExtrBuf==4) %| (SrwVisualizeWhatToExtrBuf==5)) // Phase or Re(E)
				if((SrwVisualizeWhatToExtrBuf==7) %| (SrwVisualizeWhatToExtrBuf==8)) // Phase or Re(E)
					//SrwCmpnNo=SrwVisualizeWhatToExtrBuf - 1
					SrwCmpnNo=SrwVisualizeWhatToExtrBuf - 4
				else	
					if(SrwVisualizeWhatToExtrBuf==9)
						SrwCmpnNo=7
					else
						if(SrwVisualizeWhatToExtrBuf == 2) // Flux
							if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e
								SrwCmpnNo=5
							else // Multi-e
								SrwCmpnNo=6
							endif
						else
							if(SrwVisualizeWhatToExtrBuf == 1) //Intensity
								if(SrwVisualizeSingleOrMultiEBuf==1) // Single-e
									SrwCmpnNo=1
								else // Multi-e
									SrwCmpnNo=2
								endif
							endif
						endif
					endif
				endif
				CallWfr2Int = 1
			endif
		endif
	endif
endif
if(SrwVisualizeRadTypeBuf==2) // Stokes

	if(SrwVisualizeFluxFromIntStokes == 1)
		if((SrwVisualizeWhatToExtrBuf==5) %| (SrwVisualizeWhatToExtrBuf==6))
			CallSto2PolRate = 1
			PolRateType = SrwVisualizeWhatToExtrBuf - 4
		else
			CallSto2Int = 1
		endif
	else
		if(SrwVisualizeStokesIntAndFlux==0)
			//if(SrwVisualizeWhatToExtrBuf==2) // polarization rate
			if((SrwVisualizeWhatToExtrBuf==2) %| (SrwVisualizeWhatToExtrBuf==3)) // polarization rate
				CallSto2PolRate = 1
				PolRateType = SrwVisualizeWhatToExtrBuf - 1
			else
				CallSto2Int = 1
			endif
		else
			if((SrwVisualizeWhatToExtrBuf==3) %| (SrwVisualizeWhatToExtrBuf==4)) // polarization rate
				CallSto2PolRate = 1
				PolRateType = SrwVisualizeWhatToExtrBuf - 2
			else
				CallSto2Int = 1			
			endif
		endif
		
	endif
	
endif
if(SrwVisualizeRadTypeBuf==3) // Pow Density
	CallPow2Int = 1
endif

String ComLineStr
if(SrwVisualizeRadTypeBuf==1) // Electric Field
	Variable Repr=1 // Position
	
	//if(SrwVisualizeWhatToExtrBuf != 2) // intensity, phase or Re(E) or flux
	if(CallWfr2Int == 1)
		//sprintf ComLineStr, "SrwWfr2Int(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, Repr, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		sprintf ComLineStr, "SrwWfr2Int(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, Repr, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		print ComLineStr
		SrwWfr2Int(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, Repr, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
	endif
	//if(SrwVisualizeWhatToExtrBuf==2) // polarization rate
	if(CallWfr2PolRate == 1)
		//sprintf ComLineStr, "SrwWfr2PolRate(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, Repr, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		sprintf ComLineStr, "SrwWfr2PolRateExt(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, PolRateType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		print ComLineStr
		//SrwWfr2PolRate(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, Repr, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
		SrwWfr2PolRateExt(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwCmpnNo, SrwViewPlotType, PolRateType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
	endif
endif
if(SrwVisualizeRadTypeBuf==2) // Stokes

	//if(SrwVisualizeWhatToExtrBuf != 2) // multi-e intensity
	if(CallSto2Int == 1) // multi-e intensity
		//sprintf ComLineStr, "SrwSto2Int(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		sprintf ComLineStr, "SrwSto2IntF(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, (SrwVisualizeFluxFromIntStokes + 1), SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		print ComLineStr
		//SrwSto2Int(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
		SrwSto2IntF(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, (SrwVisualizeFluxFromIntStokes + 1), SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
	endif
	//if(SrwVisualizeWhatToExtrBuf==2) // polarization rate
	if(CallSto2PolRate == 1) // polarization rate
		//sprintf ComLineStr, "SrwSto2PolRate(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		sprintf ComLineStr, "SrwSto2PolRateExt(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, (SrwVisualizeFluxFromIntStokes + 1), SrwViewPlotType, PolRateType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay
		print ComLineStr
		//SrwSto2PolRate(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, SrwViewPlotType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
		SrwSto2PolRateExt(SrwRadGenTotName, SrwSuffixExtract, SrwViewRadCmpnType, (SrwVisualizeFluxFromIntStokes + 1), SrwViewPlotType, PolRateType, SrwViewE, SrwViewX, SrwViewZ, SrwViewDisplay)
	endif
endif
//if(SrwVisualizeRadTypeBuf==3) // Pow Den
if(CallPow2Int == 1) // Pow Den
	sprintf ComLineStr, "SrwPow2Int(\"%s\",\"%s\",%g,%g,%g,%g)", SrwRadGenTotName, SrwSuffixExtract, SrwViewPlotTypeXZ, SrwViewX, SrwViewZ, SrwViewDisplay
	print ComLineStr
	SrwPow2Int(SrwRadGenTotName, SrwSuffixExtract, SrwViewPlotTypeXZ, SrwViewX, SrwViewZ, SrwViewDisplay)
endif

srwUtiSetValN("SrwVisualizeSingleOrMultiE_G", SrwVisualizeSingleOrMultiEBuf, "")
srwUtiSetValN("SrwVisualizeWhatToExtr_G", SrwVisualizeWhatToExtrBuf, "")
srwUtiSetValN("SrwVisualizeFluxCanBeExtr_G", SrwVisualizeFluxCanBeExtrBuf, "")

SrwVisualizeKillBufVars()
SrwVisualizeKillAllContr()
DoWindow/K SrwVisualizePanel

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwVisualizeHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwVisualizeDialog     ")
End
