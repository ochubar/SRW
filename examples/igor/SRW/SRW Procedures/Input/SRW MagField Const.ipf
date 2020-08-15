
//+++++++++++++++++++++++++++++++++++++++
//
//Create a Constant Magnetic field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagConstCreate(MagName,Bz,Bx)
string MagName=SrwMagConstName
Variable Bz=SrwMagConstBz
Variable Bx=SrwMagConstBx
prompt MagName,SrwPMagConstName
prompt Bz,SrwPMagConstBz
prompt Bx,SrwPMagConstBx
Silent 1						|	...
PauseUpdate

// Validation of parameters
//...

SrwMagConstName=MagName
SrwMagGenTotName=MagName+SrwMagConstType

SrwMagConstBz=Bz
SrwMagConstBx=Bx

MagName+=SrwMagConstType
Make/T/O/N=5 $MagName
//SetScale d -1E+23, 1E+23, SrwMagType_Const, $MagName

$MagName[0]=num2str(Bx)
$MagName[1]=num2str(Bz)

end

//+++++++++++++++++++++++++++++++++++++++
//
//Duplicate Constant Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagConstDupl(MagName,Name)
String MagName=SrwUndName+SrwUndType
String  Name=SrwUndName+"d"
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt Name,"Name of the Duplicated Periodic Field structure"
Silent 1						|	Duplicating the field structure  ....
PauseUpdate

SrwUndName=Name

Name += SrwUndType
duplicate/O $MagName $Name

String FullOldHarmName, MainOldHarmName, MainNewHarmName, FullNewHarmName
Variable AmOfHarm = str2num($MagName[5])
Variable i=0
do
	FullOldHarmName = $MagName[6+i]
	MainOldHarmName = FullOldHarmName[0,strlen(FullOldHarmName)-strlen(SrwUndHarmType)-2]
	MainNewHarmName = MainOldHarmName + "d"
	FullNewHarmName = MainNewHarmName + FullOldHarmName[strlen(FullOldHarmName)-strlen(SrwUndHarmType)-1] + SrwUndHarmType

	duplicate/O $FullOldHarmName $FullNewHarmName
	$Name[6+i] = FullNewHarmName
	
	i += 1
while(i<AmOfHarm)

end

//+++++++++++++++++++++++++++++++++++++++
//
//From Field Integral to Agle
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagFieldInt2Ang(FieldInt, ElecEn)
variable FieldInt=srwUtiGetValN("FieldInt", 1, "SrwUtiMagFieldInt2Ang")
variable ElecEn=SrwElecEn
prompt FieldInt,"Magnetic Field Integral [G*cm]"
prompt ElecEn,"Electron Energy [GeV]"
Silent 1						|	...
PauseUpdate

SrwElecEn=ElecEn
srwUtiSetValN("FieldInt", FieldInt, "SrwUtiMagFieldInt2Ang")

print 0.299792458e-06*FieldInt/ElecEn, "rad"
end

//+++++++++++++++++++++++++++++++++++++++
//
//Magnetic Radius
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagRad(MagField,ElecEnergy,outunit)
Variable MagField=SrwMagConstBz
Variable ElecEnergy=SrwElecEn
Variable outunit=2
prompt MagField,"Constant Magnetic Field [T]"
prompt ElecEnergy,"Electron Energy [GeV]"
prompt outunit,"Output Units",popup "mm;m;km"
Silent 1						|	...
PauseUpdate

SrwElecEn=ElecEnergy
SrwMagConstBz=MagField

String sout = ""
if(outunit == 1)
	sout = "mm"
endif
if(outunit == 2)
	sout = "m"
endif
if(outunit == 3)
	sout = "km"
endif

print srUtiMagRad(MagField, ElecEnergy, outunit), sout
end

//+++++++++++++++++++++++++++++++++++++++
//
//Critical Photon Energy / Wavelength
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagCritPhotEn(MagField,ElecEnergy,outunit)
Variable MagField=SrwMagConstBz
Variable ElecEnergy=SrwElecEn
Variable outunit=1
prompt MagField,"Constant Magnetic Field [T]"
prompt ElecEnergy,"Electron Energy [GeV]"
prompt outunit,"Output Units",popup "keV;eV;1/cm;Å;nm;µm;mm"
Silent 1						|	...
PauseUpdate

SrwElecEn=ElecEnergy
SrwMagConstBz=MagField

String sout = ""
if(outunit == 1)
	sout = "keV"
endif
if(outunit == 2)
	sout = "eV"
endif
if(outunit == 3)
	sout = "1/cm"
endif
if(outunit == 4)
	sout = "Å"
endif
if(outunit == 5)
	sout = "nm"
endif
if(outunit == 6)
	sout = "µm"
endif
if(outunit == 7)
	sout = "mm"
endif

print srUtiMagCritPhotEn(MagField, ElecEnergy, outunit), sout
end
