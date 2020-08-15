
//+++++++++++++++++++++++++++++++++++++++
//
//View SRW Structures
//
//Supports the following Types/Names:
//
//SrwElecType = "_ebm"
//SrwFieldType = "_mag"
//SrwFieldWaveType = "_fld"
//SrwUndType = "_map"
//SrwUndHarmType = "_fha"
//SrwTrjType = "_trj"
//SrwGsnBeamType = "_gsn"
//SrwSmpType = "_obs"
//SrwRadType = "_rad"
//SrwStoType = "_ras"
//SrwMomType = "_mom"
//SrwBeamlineType = "_bli"
//
//Modify necessary places (use "find"), if any new structures are supported !!!
//
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiStructTypeList()
String s=";"
String TypeList="all"+s
String/G SrwElecType, SrwFieldType, SrwFieldWaveType, SrwUndType, SrwUndHarmType, SrwTrjType
String/G SrwGsnBeamType, SrwSmpType, SrwRadType, SrwStoType, SrwMomType, SrwBeamlineType

TypeList+=srwUtiRemFirstChar(SrwElecType)+s
TypeList+=srwUtiRemFirstChar(SrwFieldType)+s
TypeList+=srwUtiRemFirstChar(SrwFieldWaveType)+s
TypeList+=srwUtiRemFirstChar(SrwUndType)+s
TypeList+=srwUtiRemFirstChar(SrwUndHarmType)+s
TypeList+=srwUtiRemFirstChar(SrwTrjType)+s
TypeList+=srwUtiRemFirstChar(SrwGsnBeamType)+s
TypeList+=srwUtiRemFirstChar(SrwSmpType)+s
TypeList+=srwUtiRemFirstChar(SrwRadType)+s
TypeList+=srwUtiRemFirstChar(SrwStoType)+s
TypeList+=srwUtiRemFirstChar(SrwMomType)+s
TypeList+=srwUtiRemFirstChar(SrwBeamlineType)+s
return TypeList
end

//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiStructTypeDescrList()
String s=";"
//String TypeList="All structures"+s
//TypeList+="Electron Beam"+s
//TypeList+="Magnetic Field"+s
//TypeList+="Magnetic Field Component"+s
//TypeList+="Periodic Magnetic Field"+s
//TypeList+="Periodic Mag. Field Harmonic"+s
//TypeList+="Trajectory"+s
//TypeList+="Gaussian Beam"+s
//TypeList+="Radiation Sampling"+s
//TypeList+="Wavefront"+s
//TypeList+="Stokes Components"+s
//TypeList+="Statistical Moments"+s
//TypeList+="Optical Component"+s

String/G SrwElecType, SrwFieldType, SrwFieldWaveType, SrwUndType, SrwUndHarmType, SrwTrjType
String/G SrwGsnBeamType, SrwSmpType, SrwRadType, SrwStoType, SrwMomType, SrwBeamlineType

String TypeList="All Types"+s
TypeList+="Electron Beam (" + SrwElecType + ")" + s
TypeList+="Magnetic Field (" + SrwFieldType + ")" + s
TypeList+="Magnetic Field Component (" + SrwFieldWaveType + ")" + s
TypeList+="Periodic Magnetic Field (" + SrwUndType + ")" + s
TypeList+="Periodic Magnetic Field Harmonic (" + SrwUndHarmType + ")" + s
TypeList+="Electron Trajectory (" + SrwTrjType + ")" + s
TypeList+="Gaussian Beam (" + SrwGsnBeamType + ")" + s
TypeList+="Radiation Sampling (" + SrwSmpType + ")" + s
TypeList+="Wavefront (" + SrwRadType + ")" + s
TypeList+="Stokes Components (" + SrwStoType + ")" + s
TypeList+="Statistical Moments (" + SrwMomType + ")" + s
TypeList+="Optical Element (" + SrwBeamlineType + ")" + s

return TypeList
end

//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiRemFirstChar(Str)
String Str
return Str[1,strlen(Str)-1]
end

//+++++++++++++++++++++++++++++++++++++++
Proc SrwUtiViewStructDialog()

SrwUtiViewStructCreateBufVars();
DoWindow/K SrwUtiViewStructPanel;
SrwUtiViewStructPanel();

End

//+++++++++++++++++++++++++++++++++++++++
Window SrwUtiViewStructPanel() : Panel

PauseUpdate; Silent 1		// building window...
//NewPanel /W=(100,50,410,200) as SrwPViewStructTitle
NewPanel /W=(100,50,428,200) as SrwPViewStructTitle

SetDrawLayer ProgBack
//DrawRect 10,13,300,110
DrawRect 10,13,316,110

SetDrawEnv fstyle= 1,fname= "default",fsize= 12
DrawText 20,33,"Structure Type and Identifier"

// Do not forget to change types here, if changed in Globals. Clean method does not work here...
PopupMenu popup0ViewStruct,value=#"srwUtiStructTypeDescrList()"   //"srwUtiStructTypeList()"; // This should match TypeList !!!
PopupMenu popup0ViewStruct,pos={25,35},size={255,21}
PopupMenu popup0ViewStruct,mode=SrwUtiViewStructListTypeNo,proc=SrwViewStructTypePopupProc

String TypeList = srwUtiStructTypeList()

//SetDrawEnv fstyle= 1,fname= "default",fsize= 12
//DrawText 140,33,SrwPViewStructType

SetDrawEnv fstyle= 1,fname= "default",fsize= 12
DrawText 20,78, "Existing Structures of this Type" //SrwPViewStructExist

Button button0ViewStruct,pos={30,120},size={70,20},proc=SrwViewStructCancelButtonProc,title=SrwPPanelCancelButton
Button button1ViewStruct,pos={130,120},size={70,20},proc=SrwViewStructOKButtonProc,title="Show"//SrwPPanelOKButton
Button button2ViewStruct,pos={231,120},size={70,20},proc=SrwViewStructHelpButtonProc,title=SrwPPanelHelpButton

String ItemName = sRWaveListItemName(TypeList, ";", SrwUtiViewStructListTypeNo)
if(cmpstr(ItemName, "") != 0)
	SrwUtiViewStructCurrentType = SrwSeparator + ItemName
endif

SrwViewStructCreatePopup()

EndMacro

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructCreatePopup()

string StrAllWaves = ""
if(cmpstr(SrwUtiViewStructCurrentType, "_all") == 0)
	//PopupMenu popup1,value=srwUtiListAllWaves()
	PopupMenu popup1,mode=1,proc=SrwViewStructNamePopupProc,value=srwUtiListAllWaves()
else
	//PopupMenu popup1,value=#"Wavelist(\"*\"+SrwUtiViewStructCurrentType,\";\",\"\")"
	PopupMenu popup1,mode=1,proc=SrwViewStructNamePopupProc,value=#"srwUtiListAllWavesValue(SrwUtiViewStructCurrentType)"
endif

PopupMenu popup1,pos={25,80},size={255,21}
//PopupMenu popup1,mode=1,proc=SrwViewStructNamePopupProc;

String AllNames
if(cmpstr(SrwUtiViewStructCurrentType, "_all") == 0)
	AllNames = Wavelist("*"+SrwElecType,";","") + Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwFieldWaveType,";","") + Wavelist("*"+SrwUndType,";","") + Wavelist("*"+SrwUndHarmType,";","");
	AllNames += Wavelist("*"+SrwTrjType,";","") + Wavelist("*"+SrwGsnBeamType,";","") + Wavelist("*"+SrwSmpType,";","") + Wavelist("*"+SrwRadType,";","") + Wavelist("*"+SrwMomType,";","") + Wavelist("*"+SrwBeamlineType,";","");
else
	AllNames = Wavelist("*"+SrwUtiViewStructCurrentType,";","");
endif
String TheName = sRWaveListItemName(AllNames, ";", 1);

//if(cmpstr(TheName, "") != 0)
SrwUtiViewStructCurrentName = TheName;
//endif

String TypeList = srwUtiStructTypeList();
String CurTypeID = srwUtiRemFirstChar(SrwUtiViewStructCurrentType);
Variable ItemNo = sRWaveListItemNo(TypeList, ";", CurTypeID);
String TypeDescrList = srwUtiStructTypeDescrList();
String TypeDescr = sRWaveListItemName(TypeDescrList, ";", ItemNo);

//SetDrawLayer/K UserBack;
//SetDrawLayer UserBack;
//DrawText 135,53,TypeDescr;

End

//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiListAllWavesValue(str)
string str

string sRes = Wavelist("*"+str,";","")
if(strlen(sRes) <= 0)
	sRes = "_none_"
endif
return sRes
end

//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiListAllWaves()
SVAR SrwElecType, SrwFieldType, SrwFieldWaveType, SrwUndType, SrwUndHarmType, SrwTrjType, SrwGsnBeamType, SrwSmpType, SrwRadType, SrwMomType, SrwBeamlineType
string strRes = Wavelist("*"+SrwElecType,";","") + Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwFieldWaveType,";","") + Wavelist("*"+SrwUndType,";","") 
strRes += Wavelist("*"+SrwUndHarmType,";","") + Wavelist("*"+SrwTrjType,";","") + Wavelist("*"+SrwGsnBeamType,";","") + Wavelist("*"+SrwSmpType,";","")
strRes += Wavelist("*"+SrwRadType,";","") + Wavelist("*"+SrwMomType,";","") + Wavelist("*"+SrwBeamlineType,";","")
return strRes
end

//+++++++++++++++++++++++++++++++++++++++
Proc SrwUtiViewStructCreateBufVars()
String/G SrwUtiViewStructCurrentType = "";
String/G SrwUtiViewStructCurrentName = "";

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwUtiViewStructKillBufVars()
KillStrings/Z SrwUtiViewStructCurrentType, SrwUtiViewStructCurrentName;

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructCancelButtonProc(ctrlName) : ButtonControl
String ctrlName;

DoWindow/K SrwUtiViewStructPanel;
SrwUtiViewStructKillBufVars();

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructOKButtonProc(ctrlName) : ButtonControl
String ctrlName;

//keep dialog alive
//DoWindow/K SrwUtiViewStructPanel;

if(cmpstr(SrwUtiViewStructCurrentName, "") != 0)
	String ExeStr;
	sprintf ExeStr, "SrwUtiViewStruct(\"%s\",%g)", SrwUtiViewStructCurrentName, 1;
	Print ExeStr;

	//SrwUtiViewStruct(SrwUtiViewStructCurrentName, 1);
	SrwUtiViewStruct(SrwUtiViewStructCurrentName, 2); //allow for modification
endif

//keep dialog alive
//SrwUtiViewStructKillBufVars();

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwUtiViewStructDialog     ")
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructTypePopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

//SrwUtiViewStructCurrentType = SrwSeparator + popStr

variable iStart = strsearch(popStr, "(", 1)
variable iEnd =strsearch(popStr, ")", iStart + 1)
SrwUtiViewStructCurrentType = popStr[iStart + 1, iEnd - 1]

KillControl popup1
SrwViewStructCreatePopup()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwViewStructNamePopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwUtiViewStructCurrentName = popStr

End

//+++++++++++++++++++++++++++++++++++++++
//
//Veiw SRW structures
//See / support structures organization here
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwUtiViewStruct(StructName, Change)
String StructName;
Variable Change;

String SrwViewType = SrwSeparator + "vew"; // Put this to globals if boss approves

variable lenStructName = strlen(StructName)
string auxStructName = ""
variable i=0
if(lenStructName > 26)
	do
		auxStructName[i] = StructName[i]
		i += 1
	while(i < 26)
else 
	auxStructName = StructName
endif

//String Names = "N" + StructName + SrwViewType;
String Names = "N" + auxStructName + SrwViewType;
sRKillAllWavesOfType(SrwViewType);

Variable Len = strlen(StructName);
String Type = StructName[Len-4] + StructName[Len-3] + StructName[Len-2] + StructName[Len-1];

if(cmpstr(Type, SrwElecType) == 0)
	Make/T/O/N=45 $Names;
	$Names[0] = "Energy [GeV]";
	$Names[1] = "Current [A]";
	$Names[2] = "<x> [m]";
	$Names[3] = "<x'> [r]";
	$Names[4] = "<z> [m]";
	$Names[5] = "<z'> [r]";
	$Names[6] = "s0 [m]";
	$Names[7] = "Hor. Emittance [nm]";
	$Names[8] = "Hor. Beta [m]";
	$Names[9] = "Hor. Alpha";
	$Names[10] = "Vert. Emittance [nm]";
	$Names[11] = "Vert. Beta [m]";
	$Names[12] = "Vert. Alpha";
	$Names[13] = "Rel. Energy Spread";
	$Names[14] = "Hor. Dispersion [m]";
	$Names[15] = "Hor. Dispersion Deriv.";
	$Names[16] = "Vert. Dispersion [m]";
	$Names[17] = "Vert. Dispersion Deriv.";
	$Names[18] = "Beam was defined from Twiss or Moments (0/1)";
	$Names[19] = "Longit. Bunch Center [m]";
	$Names[20] = "<(x-<x>)^2> [m^2]";
	$Names[21] = "<(x-<x>)(x'-<x'>)> [m]";
	$Names[22] = "<(x'-<x'>)^2> [r^2]";
	$Names[23] = "<(z-<z>)^2> [m^2]";
	$Names[24] = "<(z-<z>)(z'-<z'>)> [m]";
	$Names[25] = "<(z'-<z'>)^2> [r^2]";
	$Names[26] = "<(x-<x>)(z-<z>)> [m^2]";
	$Names[27] = "<(x'-<x'>)(z-<z>)> [m]";
	$Names[28] = "<(x-<x>)(z'-<z'>)> [m]";
	$Names[29] = "<(x'-<x'>)(z'-<z'>)> [r^2]";
	$Names[33] = "<(s-<s>)^2> [m^2]";
	$Names[34] = "<(s-<s>)(E-<E>)>/<E> [m]";
	$Names[35] = "<(x-<x>)(E-<E>)>/<E> [m]";
	$Names[36] = "<(x'-<x'>)(E-<E>)>/<E> [r^2]";
	$Names[37] = "<(z-<z>)(E-<E>)>/<E> [m]";
	$Names[38] = "<(z'-<z'>)(E-<E>)>/<E> [r^2]";
	$Names[39] = "<(x-<x>)(s-<s>)> [m^2]";
	$Names[40] = "<(x'-<x'>)(s-<s>)> [m]";
	$Names[41] = "<(z-<z>)(s-<s>)> [m^2]";
	$Names[42] = "<(z'-<z'>)(s-<s>)> [m]";
	$Names[43] = "Number of Electrons in Bunch";
endif
if(cmpstr(Type, SrwFieldType) == 0)
	Make/T/O/N=8 $Names;
	$Names[0] = "Hor. Field Wave";
	$Names[1] = "Vert. Field Wave";
	$Names[2] = "Meth. of Integration";
	$Names[3] = "Rel. Prec./Step Size [m]";
	$Names[4] = "Init. Integ. Coord. [m]";
	$Names[5] = "Fin. Integ. Coord. [m]";
	$Names[6] = "";
	$Names[7] = "Points to Save";
endif
if(cmpstr(Type, SrwFieldWaveType) == 0)
	// Do nothing
endif
if(cmpstr(Type, SrwUndType) == 0)
	Make/T/O/N=(6 + str2num($StructName[5])) $Names;
	$Names = "Harmonic";
	$Names[0] = "Period [m]";
	$Names[1] = "Length [m]";
	$Names[2] = "Type (1/2/3/4)";
	$Names[3] = "";
	$Names[4] = "Fundamental [keV/GeV^2]";
	$Names[5] = "Number of Harmonics";
endif
if(cmpstr(Type, SrwUndHarmType) == 0)
	Make/T/O/N=4 $Names;
	$Names[0] = "Harmonic";
	$Names[1] = "Plane Vert./Hor. (1/2)";
	$Names[2] = "K";
	$Names[3] = "Phase [r]";
endif
if(cmpstr(Type, SrwTrjType) == 0)
	Make/T/O/N=6 $Names;
	$Names[0] = "Hor. Pos. vs Long. Pos. wave";
	$Names[1] = "Vert. Pos. vs Long. Pos. wave";
	$Names[2] = "";
	$Names[3] = "";
	$Names[4] = "Energy [GeV]";
	$Names[5] = "Current [A]";
endif
if(cmpstr(Type, SrwGsnBeamType) == 0)
	Make/T/O/N=7 $Names;
	$Names[0] = "Elec. Beam struct.";
	$Names[1] = "Hor. Waist RMS Size [m]";
	$Names[2] = "Vert. Waist RMS Size [m]";
	$Names[3] = "Hor. Mode Order";
	$Names[4] = "Vert. Mode Order";
	$Names[5] = "Photons/(0.1%BW)";
	$Names[6] = "Polarization";
endif
if(cmpstr(Type, SrwSmpType) == 0)
	Make/T/O/N=14 $Names;
	$Names[0] = "";
	$Names[1] = "";
	$Names[2] = "";
	$Names[3] = "";
	$Names[4] = "Longitud. Pos. [m]";
	$Names[5] = "Init. Photon Energy [eV]";
	$Names[6] = "Fin. Photon Energy [eV]";
	$Names[7] = "Number of Energy Points";
	$Names[8] = "Init. Hor. Pos. [m]";
	$Names[9] = "Fin. Hor. Pos. [m]";
	$Names[10] = "Number of Hor. Points";
	$Names[11] = "Init. Vert. Pos. [m]";
	$Names[12] = "Fin. Vert. Pos. [m]";
	$Names[13] = "Number of Vert. Points";
endif
if(cmpstr(Type, SrwRadType) == 0)
	Make/T/O/N=22 $Names;
	$Names[0] = "Hor. Electric Field";
	$Names[1] = "Vert. Electric Field";
	$Names[2] = "Coor./Ang. Repr. (0/1)";
	$Names[3] = "Step of Phot. Energy [eV]";
	$Names[4] = "Init. Phot. Energy [eV]";
	$Names[5] = "Step of Hor. Pos. [m]";
	$Names[6] = "Init. Hor. Pos. [m]";
	$Names[7] = "Step of Vert. Pos. [m]";
	$Names[8] = "Init. Vert. Pos. [m]";
	$Names[9] = "Quality of Propag.";
	$Names[10] = "Freq./Time Repr. (0/1)";
	$Names[11] = "";
	$Names[12] = "Elec. Beam struct.";
	$Names[13] = "Moments Transf. Matr.";
	$Names[14] = "Linearity of Transf.";
	$Names[15] = "Rad. Moments Hor. Pol.";
	$Names[16] = "Rad. Moments Vert. Pol.";
	$Names[17] = "Moments Pred./Comp. (0/1)";
	$Names[18] = "Auxiliary Wavefront Data";
	$Names[19] = "";
	$Names[20] = "";
	$Names[21] = "";
endif
if(cmpstr(Type, SrwStoType) == 0)
	// Do nothing
endif
if(cmpstr(Type, SrwMomType) == 0)
	Make/T/O/N=11 $Names;
	$Names[0] = "Total Phot./s/(0.1%BW)";
	$Names[1] = "<x> [m]";
	$Names[2] = "<x'> [r]";
	$Names[3] = "<z> [m]";
	$Names[4] = "<z'> [r]";
	$Names[5] = "<x^2> [m^2]";
	$Names[6] = "<xx'> [m]";
	$Names[7] = "<x'^2> [r^2]";
	$Names[8] = "<z^2> [m^2]";
	$Names[9] = "<zz'> [m]";
	$Names[10] = "<z'^2> [r^2]";
endif
if(cmpstr(Type, SrwBeamlineType) == 0)
	String ElemIDStr = $StructName[0];
	
	if(cmpstr(ElemIDStr, SrwBliDriftType) == 0)
		Make/T/O/N=2 $Names;
		$Names[0] = "Element Type";
		$Names[1] = "Length [m]";
	endif
	if(cmpstr(ElemIDStr, SrwBliThinLensType) == 0)
		Make/T/O/N=5 $Names;
		$Names[0] = "Element Type";
		$Names[1] = "Hor. Focal Length [m]";
		$Names[2] = "Vert. Focal Length [m]";
		$Names[3] = "Hor. Center Pos. [m]";
		$Names[4] = "Vert. Center Pos. [m]";
	endif
	if(cmpstr(ElemIDStr, SrwBliRectApertType) == 0)
		Make/T/O/N=5 $Names;
		$Names[0] = "Element Type";
		$Names[1] = "Hor. Size [m]";
		$Names[2] = "Vert. Size [m]";
		$Names[3] = "Hor. Center Pos. [m]";
		$Names[4] = "Vert. Center Pos. [m]";
	endif
	if(cmpstr(ElemIDStr, SrwBliCircApertType) == 0)
		Make/T/O/N=4 $Names;
		$Names[0] = "Element Type";
		$Names[1] = "Diameter [m]";
		$Names[2] = "Hor. Center Pos. [m]";
		$Names[3] = "Vert. Center Pos. [m]";
	endif
	if(cmpstr(ElemIDStr, SrwBliSpherMirrorType) == 0)
		Make/T/O/N=8 $Names;
		$Names[0] = "Element Type";
		$Names[1] = "Radius of Curvature [m]";
		$Names[2] = "Hor. Size [m]";
		$Names[3] = "Vert. Size [m]";
		$Names[4] = "Rotation Plane";
		$Names[5] = "Rotation Angle [r]";
		$Names[6] = "Hor. Center Pos. [m]";
		$Names[7] = "Vert. Center Pos. [m]";
	endif
	if(cmpstr(ElemIDStr, SrwBliContType) == 0)
		Make/T/O/N=(DimSize($StructName,0)) $Names;
		$Names = "Member Element";
		$Names[0] = "Element Type";
	endif
	// Process more optical Elements here when they appear
endif

String pStructToView = StructName + SrwViewType;

if(Change == 1) // Do not change wave
	Duplicate/O $StructName $pStructToView;
else
	pStructToView = StructName;
endif

Variable StructExists = WaveExists($pStructToView);
Variable NamesExist = WaveExists($Names);

if(StructExists == 1)
	if(NamesExist == 1)
		Edit $Names, $pStructToView;
	else
		Edit $pStructToView.id;
	endif
endif

End

//+++++++++++++++++++++++++++++++++++++++
