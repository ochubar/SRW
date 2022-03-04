@echo off
REM -- First make map file from Microsoft Visual C++ generated resource.h
echo // MAKEHELP.BAT generated Help Map file.  Used by SETUP.HPJ. >"hlp\Setup.hm"
echo. >>"hlp\Setup.hm"
echo // Commands (ID_* and IDM_*) >>"hlp\Setup.hm"
makehm ID_,HID_,0x10000 IDM_,HIDM_,0x10000 resource.h >>"hlp\Setup.hm"
echo. >>"hlp\Setup.hm"
echo // Prompts (IDP_*) >>"hlp\Setup.hm"
makehm IDP_,HIDP_,0x30000 resource.h >>"hlp\Setup.hm"
echo. >>"hlp\Setup.hm"
echo // Resources (IDR_*) >>"hlp\Setup.hm"
makehm IDR_,HIDR_,0x20000 resource.h >>"hlp\Setup.hm"
echo. >>"hlp\Setup.hm"
echo // Dialogs (IDD_*) >>"hlp\Setup.hm"
makehm IDD_,HIDD_,0x20000 resource.h >>"hlp\Setup.hm"
echo. >>"hlp\Setup.hm"
echo // Frame Controls (IDW_*) >>"hlp\Setup.hm"
makehm IDW_,HIDW_,0x50000 resource.h >>"hlp\Setup.hm"
REM -- Make help for Project SETUP


echo Building Win32 Help files
start /wait hcw /C /E /M "hlp\Setup.hpj"
if errorlevel 1 goto :Error
if not exist "hlp\Setup.hlp" goto :Error
if not exist "hlp\Setup.cnt" goto :Error
echo.
if exist Debug\nul copy "hlp\Setup.hlp" Debug
if exist Debug\nul copy "hlp\Setup.cnt" Debug
if exist Release\nul copy "hlp\Setup.hlp" Release
if exist Release\nul copy "hlp\Setup.cnt" Release
echo.
goto :done

:Error
echo hlp\Setup.hpj(1) : error: Problem encountered creating help file

:done
echo.
