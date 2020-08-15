
//+++++++++++++++++++++++++++++++++++++++
//
// Report Generation Utilities
//
//+++++++++++++++++++++++++++++++++++++++
//Add frame and grid to current graph
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiGraphAddFrameAndGrid()
if(strlen(WinName(0,1)) > 0)
	ModifyGraph grid(left)=1,tick(left)=2,mirror(left)=1
	ModifyGraph grid=1,tick=2,mirror=1
	ModifyGraph notation(left)=1
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Add axes labels to current graph
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiGraphAddLabels(VertLabel, HorLabel)
string HorLabel=srwUtiGetValS("HorLabel", "", "SrwUtiGraphAddLabels")
string VertLabel=srwUtiGetValS("VertLabel", "", "SrwUtiGraphAddLabels")
prompt HorLabel,"Horizontal (bottom) label"
prompt VertLabel,"Vertical (left) label"
srwUtiSetValS("HorLabel", HorLabel, "SrwUtiGraphAddLabels")
srwUtiSetValS("HorLabel", HorLabel, "SrwUtiGraphAddLabels")

Label bottom HorLabel
Label left VertLabel
end

//+++++++++++++++++++++++++++++++++++++++
//Resize and position current graph window
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiGraphWindResize(x0, y0, xSize, ySize, leftAxOffset, bottomAxOffset)
variable x0=srwUtiGetValN("x0", 10, "SrwUtiGraphWindResize")
variable y0=srwUtiGetValN("y0", 10, "SrwUtiGraphWindResize")
variable xSize=srwUtiGetValN("xSize", 200, "SrwUtiGraphWindResize")
variable ySize=srwUtiGetValN("ySize", 150, "SrwUtiGraphWindResize")
variable leftAxOffset=srwUtiGetValN("leftAxOffset", 0, "SrwUtiGraphWindResize")
variable bottomAxOffset=srwUtiGetValN("bottomAxOffset", 0, "SrwUtiGraphWindResize")
prompt x0,"Horizontal position of top-left corner"
prompt y0,"Vertical position of top-left corner"
prompt xSize,"Horizontal window size"
prompt ySize,"Vertical window size"
prompt leftAxOffset,"Left axis offset (-10..10)"
prompt bottomAxOffset,"Bottom axis offset (-10..10)"
srwUtiSetValN("x0", x0, "SrwUtiGraphWindResize")
srwUtiSetValN("y0", y0, "SrwUtiGraphWindResize")
srwUtiSetValN("xSize", xSize, "SrwUtiGraphWindResize")
srwUtiSetValN("ySize", ySize, "SrwUtiGraphWindResize")
srwUtiSetValN("leftAxOffset", leftAxOffset, "SrwUtiGraphWindResize")
srwUtiSetValN("bottomAxOffset", bottomAxOffset, "SrwUtiGraphWindResize")

variable xe = x0 + xSize
variable ye = y0 + ySize
MoveWindow x0, y0, xe, ye
ModifyGraph axOffset(left)=leftAxOffset
ModifyGraph axOffset(bottom)=bottomAxOffset
end

//+++++++++++++++++++++++++++++++++++++++
//Save graph to a graphics file (PNG only)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiGraphSave(FolderPath, FileName)
string FolderPath=srwUtiGetValS("FolderPath", "", "SrwUtiGraphSave")
string FileName=srwUtiGetValS("FileName", "", "SrwUtiGraphSave")
prompt FolderPath,"Path to the folder where to save the file"
prompt FileName,"Graphics file name without extension"
srwUtiSetValS("FolderPath", FolderPath, "SrwUtiGraphSave")
srwUtiSetValS("FileName", FileName, "SrwUtiGraphSave")

if(strlen(FolderPath) > 0)
	NewPath/O/Q pathName, FolderPath
	if(strlen(FileName) > 0)
		string TotFileName = FileName + ".png"
		SavePICT/O/E=-5/B=144/P=pathName as TotFileName
	else
		SavePICT/O/E=-5/B=144/P=pathName
	endif
else
	if(strlen(FileName) > 0)
		string TotFileName = FileName + ".png"
		SavePICT/O/E=-5/B=144 as TotFileName
	else
		SavePICT/O/E=-5/B=144
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Save graph to a graphics file (several formats supported)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiGraphSaveFile(FolderPath, FileName, GraphFormatNum)
string FolderPath=srwUtiGetValS("FolderPath", "", "SrwUtiGraphSave")
string FileName=srwUtiGetValS("FileName", "", "SrwUtiGraphSave")
variable GraphFormatNum=srwUtiGetValN("GraphFormatNum", 1, "SrwUtiGraphSave")
prompt FolderPath,"Path to the folder where to save the file"
prompt FileName,"Graphics file name without extension"
prompt GraphFormatNum,"Graphic file format",popup "PNG;JPEG;TIFF;PDF;EPS"
srwUtiSetValS("FolderPath", FolderPath, "SrwUtiGraphSave")
srwUtiSetValS("FileName", FileName, "SrwUtiGraphSave")
srwUtiSetValN("GraphFormatNum", GraphFormatNum, "SrwUtiGraphSave")

make/O wIgorGraphFormID = {-5,-6,-7,-8,-3} //"PNG;JPEG;TIFF;PDF;EPS"
make/O/T wGraphFileExt = {"png","jpg","tiff","pdf","eps"}

variable igorGraphFormID = wIgorGraphFormID[GraphFormatNum - 1]
string fileExt = wGraphFileExt[GraphFormatNum - 1]

if(strlen(FolderPath) > 0)
	NewPath/O/Q pathName, FolderPath
	if(strlen(FileName) > 0)
		string TotFileName = FileName + "." + fileExt
		//SavePICT/O/E=-5/B=144/P=pathName as TotFileName
		SavePICT/O/E=(igorGraphFormID)/B=144/P=pathName as TotFileName
	else
		//SavePICT/O/E=-5/B=144/P=pathName
		SavePICT/O/E=(igorGraphFormID)/B=144/P=pathName
	endif
else
	if(strlen(FileName) > 0)
		string TotFileName = FileName + "." + fileExt
		//SavePICT/O/E=-5/B=144 as TotFileName
		SavePICT/O/E=(igorGraphFormID)/B=144 as TotFileName
	else
		//SavePICT/O/E=-5/B=144
		SavePICT/O/E=(igorGraphFormID)/B=144
	endif
endif

killwaves/Z wIgorGraphFormID, wGraphFileExt
end