#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//Soleil IDs linked to straight sections
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiSrcSelect(src)
variable src = srwUtiGetValN("src", 2, "SrwUtiSrcSelect")
prompt src, "Select Synchrotron Radiation Source", popup SrwUtiSrcListAvail(";")
Silent 1						|	...
PauseUpdate
srwUtiSetValN("src", src, "SrwUtiSrcSelect")

variable/G SrwUtiSrcSelected = src
string strSrcSel = StringFromList(src - 1, SrwUtiSrcListAvail(";"))

//Define text wave with ID names and their corresponding file names
string nmGWID = "gwSrwUtiSrcIDs", nmGWS = "gwSrwUtiSrcStrSect", nmGWMP = "gwSrwUtiSrcMainParam", nmGWBM = "gwSrwUtiSrcBM"
make/O/T/N=(1,3) $nmGWID
make/O/T/N=(1,3) $nmGWS
make/O/T/N=3 $nmGWMP
make/O/T/N=(2,5) $nmGWBM

if(cmpstr(strSrcSel, "SOLEIL") == 0)
	redimension/N=(17, 4) $nmGWID
	$nmGWID[0][0] = "U20"; 							$nmGWID[0][1] = ""; 									$nmGWID[0][2] = "Bz_U20.dta"; 						$nmGWID[0][3] = "0.2"
	$nmGWID[1][0] = "HU45 Linear Horizontal Polar."; 	$nmGWID[1][1] = ""; 									$nmGWID[1][2] = "Bz_HU45_LinHor.dta"; 				$nmGWID[1][3] = "0.9225"
	$nmGWID[2][0] = "HU45 Helical Mode"; 			$nmGWID[2][1] = "Bx_HU45_Hel.dta";					$nmGWID[2][2] = "Bz_HU45_Hel.dta"; 					$nmGWID[2][3] = "0.9225"
	$nmGWID[3][0] = "HU45 Linear Vertical Polar."; 		$nmGWID[3][1] = "Bx_HU45_LinVert.dta"; 				$nmGWID[3][2] = "Bz_HU45_LinVert.dta"; 				$nmGWID[3][3] = "0.9225"
	$nmGWID[4][0] = "HU60 Linear Horizontal Polar."; 	$nmGWID[4][1] = ""; 									$nmGWID[4][2] = "Bz_HU60_15_LinHor.dta"; 			$nmGWID[4][3] = "1.05"
	$nmGWID[5][0] = "HU60 Helical Mode"; 			$nmGWID[5][1] = "Bx_HU60_15_Hel.dta"; 				$nmGWID[5][2] = "Bz_HU60_15_Hel.dta"; 				$nmGWID[5][3] = "1.05"
	$nmGWID[6][0] = "HU60 Linear Vertical Polar."; 		$nmGWID[6][1] = "Bx_HU60_15_LinVert.dta"; 			$nmGWID[6][2] = "Bz_HU60_15_LinVert.dta"; 			$nmGWID[6][3] = "1.05"
	$nmGWID[7][0] = "HU80 Linear Horizontal Polar."; 	$nmGWID[7][1] = ""; 									$nmGWID[7][2] = "Bz_HU80_ELETTRA_15_LinHor.dta"; 	$nmGWID[7][3] = "1.08"
	$nmGWID[8][0] = "HU80 Helical Mode"; 			$nmGWID[8][1] = "Bx_HU80_ELETTRA_15_Hel.dta"; 	$nmGWID[8][2] = "Bz_HU80_ELETTRA_15_Hel.dta"; 	$nmGWID[8][3] = "1.08"
	$nmGWID[9][0] = "HU80 Linear Vertical Polar."; 		$nmGWID[9][1] = "Bx_HU80_ELETTRA_15_LinVert.dta"; 	$nmGWID[9][2] = "Bz_HU80_ELETTRA_15_LinVert.dta"; 	$nmGWID[9][3] = "1.08"
	$nmGWID[10][0] = "HU256 Linear Horizontal Polar."; 	$nmGWID[10][1] = "Bx_HU256_0A2650A.dta"; 			$nmGWID[10][2] = "Bz_HU256_0A2650A.dta"; 			$nmGWID[10][3] = "1.28"
	$nmGWID[11][0] = "HU256 Helical Mode"; 			$nmGWID[11][1] = "Bx_HU256_Hel.dta"; 				$nmGWID[11][2] = "Bz_HU256_Hel.dta"; 				$nmGWID[11][3] = "1.28"
	$nmGWID[12][0] = "HU256 Linear Vertical Polar."; 	$nmGWID[12][1] = "Bx_HU256_Hel.dta"; 				$nmGWID[12][2] = ""; 								$nmGWID[12][3] = "1.28"
	$nmGWID[13][0] = "HU640 Linear Horizontal Polar."; 	$nmGWID[13][1] = ""; 								$nmGWID[13][2] = "Bz_HU640_Plan.dta"; 				$nmGWID[13][3] = "6.4"
	$nmGWID[14][0] = "HU640 Helical Mode"; 			$nmGWID[14][1] = "Bx_HU640_Plan.dta"; 				$nmGWID[14][2] = "Bz_HU640_Plan.dta"; 				$nmGWID[14][3] = "6.4"
	$nmGWID[15][0] = "HU640 Linear Vertical Polar."; 	$nmGWID[15][1] = "Bx_HU640_Plan.dta"; 				$nmGWID[15][2] = ""; 								$nmGWID[15][3] = "6.4"
	$nmGWID[16][0] = "HU640 Linear Tilted Polar."; 		$nmGWID[16][1] = "Bx_HU640_Plan.dta"; 				$nmGWID[16][2] = "Bz_HU640_Tilt.dta"; 				$nmGWID[16][3] = "6.4"
	
	redimension/N=(4, 3) $nmGWS
	$nmGWS[0][0] = "Short";						$nmGWS[0][1] = "7.722";								$nmGWS[0][2] = "0.001016,3.73,0.037,17.78,1.75,0,0,0.28,0"; //thick beam params
	$nmGWS[1][0] = "Medium";					$nmGWS[1][1] = "12.405";							$nmGWS[1][2] = "0.001016,3.73,0.037,4,1.77,0,0,0.13,0";
	$nmGWS[2][0] = "Long";						$nmGWS[2][1] = "18.54";								$nmGWS[2][2] = "0.001016,3.73,0.037,10.09,8.01,0,0,0.2,0";
	$nmGWS[3][0] = "Very Short";					$nmGWS[3][1] = "4.4527";							$nmGWS[3][2] = "0.001016,3.73,0.037,17.78,1.75,0,0,0.28,0"
	
	redimension/N=(2,5) $nmGWBM
	$nmGWBM[0][0] = "No";						$nmGWBM[0][1] = "";									$nmGWBM[0][2] = "";									$nmGWBM[0][3] = "";		$nmGWBM[0][4] = ""
	$nmGWBM[1][0] = "Yes";						$nmGWBM[1][1] = "BzSolDipole.dta";					$nmGWBM[1][2] = "BzSolDipole.dta";					$nmGWBM[1][3] = "5";	$nmGWBM[1][4] = "0.52537"

	redimension/N=2 $nmGWMP
	$nmGWMP[0] = "SOLEIL"	;					$nmGWMP[1] = "2.75,0.5,0,0,0,0,0";	//$nmGWMP[2] = "BzSolDipole.dta"
endif
if((cmpstr(strSrcSel, "NSLS-II") == 0) %| (cmpstr(strSrcSel, "NSLSII") == 0) %| (cmpstr(strSrcSel, "NSLS-2") == 0) %| (cmpstr(strSrcSel, "NSLS2") == 0))
	redimension/N=(38, 4) $nmGWID
	$nmGWID[0][0] = "IVU20 (5 mm gap)";									$nmGWID[0][1] = ""; 										$nmGWID[0][2] = "IVU20G5CBZ.dat"; 						$nmGWID[0][3] = "0.206680445363"
	$nmGWID[1][0] = "IVU21 (5.5 mm gap)";								$nmGWID[1][1] = ""; 										$nmGWID[1][2] = "IVU21G5_5CBZ.dat";						$nmGWID[1][3] = "0.21"
	$nmGWID[2][0] = "IVU22 (6.95 mm gap)";								$nmGWID[2][1] = ""; 										$nmGWID[2][2] = "IVU22G6_95CBZ.dat"; 					$nmGWID[2][3] = "0.22"
	//$nmGWID[3][0] = "EPU49 (2 m) Linear Horizontal Polar. Mode";			$nmGWID[3][1] = ""; 										$nmGWID[3][2] = "EPU49_BZ_LH30x30_40perC.dat"; 		$nmGWID[3][3] = "0.49"
	$nmGWID[3][0] = "EPU49 (2 m) Linear Horizontal Polar. Mode";			$nmGWID[3][1] = ""; 										$nmGWID[3][2] = "EPU49G11_5LHCBZ.dat"; 				$nmGWID[3][3] = "0.0492"
	$nmGWID[4][0] = "EPU49 (2 m) Helical Mode";							$nmGWID[4][1] = "EPU49_BX_HE30x30_40perC.dat"; 		$nmGWID[4][2] = "EPU49_BZ_HE30x30_40perC.dat"; 		$nmGWID[4][3] = "0.49"
	$nmGWID[5][0] = "EPU49 (2 m) Linear Vertical Polar. Mode";				$nmGWID[5][1] = "EPU49_BX_LV30x30_40perC.dat"; 		$nmGWID[5][2] = "EPU49_BZ_LV30x30_40perC.dat"; 		$nmGWID[5][3] = "0.49"
	$nmGWID[6][0] = "EPU49 (2 m) Linear Tilted (45 deg.) Polar. Mode";		$nmGWID[6][1] = "EPU49_BX_LT45_30x30_40perC.dat"; 		$nmGWID[6][2] = "EPU49_BZ_LT45_30x30_40perC.dat"; 		$nmGWID[6][3] = "0.49"
	$nmGWID[7][0] = "2 x EPU49 Inline Linear Horizontal Polar. Mode";		$nmGWID[7][1] = ""; 										$nmGWID[7][2] = "InlineEPU49_BZ_LH30x30_40per.dat"; 		$nmGWID[7][3] = "0.13" //"0.5200520052"
	$nmGWID[8][0] = "2 x EPU49 Inline Helical Mode";						$nmGWID[8][1] = "InlineEPU49_BX_HE30x30_40per.dat";		$nmGWID[8][2] = "InlineEPU49_BZ_HE30x30_40per.dat"; 	$nmGWID[8][3] = "0.13" //"0.5200520052"
	$nmGWID[9][0] = "2 x EPU49 Inline Linear Vertical Polar. Mode";			$nmGWID[9][1] = "InlineEPU49_BX_LV30x30_40per.dat"; 		$nmGWID[9][2] = "InlineEPU49_BZ_LV30x30_40per.dat"; 		$nmGWID[9][3] = "0.1238" //"0.24761237585"
	$nmGWID[10][0] = "2 x EPU49 Inline Linear Tilted (45 deg.) Polar. Mode";	$nmGWID[10][1] = "InlineEPU49_BX_LT45_30x30_40per.dat"; 	$nmGWID[10][2] = "InlineEPU49_BZ_LT45_30x30_40per.dat"; 	$nmGWID[10][3] = "0.1238" //"0.24761238"
	$nmGWID[11][0] = "2 x EPU49 Canted Linear Horizontal Polar. Mode";		$nmGWID[11][1] = ""; 									$nmGWID[11][2] = "CantEPU49_BZ_LH30x30_40per.dat";		$nmGWID[11][3] = "0.56"
	$nmGWID[12][0] = "2 x EPU49 Canted Helical Mode";					$nmGWID[12][1] = "CantEPU49_BX_HE30x30_40per.dat";	$nmGWID[12][2] = "CantEPU49_BZ_HE30x30_40per.dat";	$nmGWID[12][3] = "0.29"
	$nmGWID[13][0] = "2 x EPU49 Canted Linear Vertical Polar. Mode";		$nmGWID[13][1] = "CantEPU49_BX_LV30x30_40per.dat";		$nmGWID[13][2] = "CantEPU49_BZ_LV30x30_40per.dat"; 	$nmGWID[13][3] = "0.5800580058"
	$nmGWID[14][0] = "2 x EPU49 Canted Linear Tilted (45 deg.) Polar. Mode";	$nmGWID[14][1] = "CantEPU49_BX_LT45_30x30_40per.dat";	$nmGWID[14][2] = "CantEPU49_BZ_LT45_30x30_40per.dat";	$nmGWID[14][3] = "0.5800580058"
	$nmGWID[15][0] = "EPU49 (4 m) Linear Horizontal Polar. Mode";			$nmGWID[15][1] = "";										$nmGWID[15][2] = "EPU49_BZ_LH30x30_81perC.dat"; 		$nmGWID[15][3] = "0.49"
	$nmGWID[16][0] = "EPU49 (4 m) Helical Mode";							$nmGWID[16][1] = "EPU49_BX_HE30x30_81perC.dat"; 		$nmGWID[16][2] = "EPU49_BZ_HE30x30_81perC.dat"; 		$nmGWID[16][3] = "0.49"
	$nmGWID[17][0] = "EPU49 (4 m) Linear Vertical Polar. Mode";				$nmGWID[17][1] = "EPU49_BX_LV30x30_81perC.dat"; 		$nmGWID[17][2] = "EPU49_BZ_LV30x30_81perC.dat"; 		$nmGWID[17][3] = "0.49"
	$nmGWID[18][0] = "EPU49 (4 m) Linear Tilted (45 deg.) Polar. Mode";		$nmGWID[18][1] = "EPU49_BX_LT45_30x30_81perC.dat"; 	$nmGWID[18][2] = "EPU49_BZ_LT45_30x30_81perC.dat"; 	$nmGWID[18][3] = "0.49"
	$nmGWID[19][0] = "EPU105 (2.8 m, 16 mm gap) Lin. Horizontal Polar. Mode";	$nmGWID[19][1] = "";									$nmGWID[19][2] = "EPU105L2_8G16_BZ_LHc.dat"; 			$nmGWID[19][3] = "0.1052"
	$nmGWID[20][0] = "EPU105 (2.8 m, 16 mm gap) Helical Mode";			$nmGWID[20][1] = "EPU105L2_8G16_BX_HEc.dat"; 			$nmGWID[20][2] = "EPU105L2_8G16_BZ_HEc.dat"; 			$nmGWID[20][3] = "0.1052"
	$nmGWID[21][0] = "EPU105 (2.8 m, 22 mm gap) Helical Mode";			$nmGWID[21][1] = "EPU105L2_8G22_BX_HEc.dat"; 			$nmGWID[21][2] = "EPU105L2_8G22_BZ_HEc.dat"; 			$nmGWID[21][3] = "0.1052"
	$nmGWID[22][0] = "EPU105 (2.8 m, 16 mm gap) Lin. Vertical Polar. Mode";	$nmGWID[22][1] = "EPU105L2_8G16_BX_LVc.dat"; 			$nmGWID[22][2] = "EPU105L2_8G16_BZ_LVc.dat"; 			$nmGWID[22][3] = "0.1052"	
	$nmGWID[23][0] = "EPU105 (2.8 m, 19 mm gap) Lin. Vertical Polar. Mode";	$nmGWID[23][1] = "EPU105L2_8G19_BX_LVc.dat"; 			$nmGWID[23][2] = "EPU105L2_8G19_BZ_LVc.dat"; 			$nmGWID[23][3] = "0.1052"	
	$nmGWID[24][0] = "EPU105 (2.8 m, 16 mm gap) Lin. Tilted (45 deg.) Polar. Mode";	$nmGWID[24][1] = "EPU105L2_8G16_BX_LTc.dat"; 	$nmGWID[24][2] = "EPU105L2_8G16_BZ_LTc.dat"; 			$nmGWID[24][3] = "0.1052"
	$nmGWID[25][0] = "EPU57 (1.4 m, 16 mm gap) Lin. Horizontal Polar. Mode";$nmGWID[25][1] = "";									$nmGWID[25][2] = "EPU57L1_4G16_BZ_LHc.dat"; 			$nmGWID[25][3] = "0.0572"
	$nmGWID[26][0] = "EPU57 (1.4 m, 16 mm gap) Helical Mode";			$nmGWID[26][1] = "EPU57L1_4G16_BX_HEc.dat"; 			$nmGWID[26][2] = "EPU57L1_4G16_BZ_HEc.dat"; 			$nmGWID[26][3] = "0.0572"
	$nmGWID[27][0] = "EPU57 (1.4 m, 16 mm gap) Lin. Vertical Polar. Mode";	$nmGWID[27][1] = "EPU57L1_4G16_BX_LVc.dat"; 			$nmGWID[27][2] = "EPU57L1_4G16_BZ_LVc.dat"; 			$nmGWID[27][3] = "0.0572"
	$nmGWID[28][0] = "EPU57 (1.4 m, 16 mm gap) Lin. Tilted (45 deg.) Polar. Mode";	$nmGWID[28][1] = "EPU57L1_4G16_BX_LTc.dat"; 	$nmGWID[28][2] = "EPU57L1_4G16_BZ_LTc.dat"; 			$nmGWID[28][3] = "0.0572"
	$nmGWID[29][0] = "EPU57 (7 m, 16 mm gap) Lin. Horizontal Polar. Mode"; $nmGWID[29][1] = "";										$nmGWID[29][2] = "EPU57L7G16_BZ_LHc.dat"; 				$nmGWID[29][3] = "0.0572"
	$nmGWID[30][0] = "EPU57 (7 m, 16 mm gap) Helical Mode";				$nmGWID[30][1] = "EPU57L7G16_BX_HEc.dat"; 			$nmGWID[30][2] = "EPU57L7G16_BZ_HEc.dat"; 			$nmGWID[30][3] = "0.0572"	
	$nmGWID[31][0] = "EPU57 (7 m, 16 mm gap) Lin. Vertical Polar. Mode";	$nmGWID[31][1] = "EPU57L7G16_BX_LVc.dat"; 				$nmGWID[31][2] = "EPU57L7G16_BZ_LVc.dat"; 				$nmGWID[31][3] = "0.0572"
	$nmGWID[32][0] = "EPU57 (7 m, 16 mm gap) Lin. Tilted (45 deg.) Polar. Mode";	$nmGWID[32][1] = "EPU57L7G16_BX_LTc.dat"; 			$nmGWID[32][2] = "EPU57L7G16_BZ_LTc.dat"; 				$nmGWID[32][3] = "0.0572"
	$nmGWID[33][0] = "DW100";											$nmGWID[33][1] = ""; 									$nmGWID[33][2] = "DW100CBZ.dat"; 						$nmGWID[33][3] = "1.49253731343284" //"1"
	$nmGWID[34][0] = "2 x DW100 Inline";									$nmGWID[34][1] = ""; 									$nmGWID[34][2] = "Inline2DW100BZ.dat"; 					$nmGWID[34][3] = "0.4" //"1"
	$nmGWID[35][0] = "DW90";											$nmGWID[35][1] = ""; 									$nmGWID[35][2] = "DW90CBZ.dat"; 						$nmGWID[35][3] = "1"
	$nmGWID[36][0] = "2 x DW90 Inline";									$nmGWID[36][1] = ""; 									$nmGWID[36][2] = "Inline2DW90BZ.dat"; 					$nmGWID[36][3] = "0.191754796254"
	$nmGWID[37][0] = "TPW Isolated";										$nmGWID[37][1] = ""; 									$nmGWID[37][2] = "TPWCBZ.dat"; 							$nmGWID[37][3] = "0.108"
	//$nmGWID[18][0] = "TPW with Bending Magnet Edges";					$nmGWID[18][1] = ""; 									$nmGWID[18][2] = "TPWwithBMBZ.dat"; 					$nmGWID[18][3] = "0.1800045"

	redimension/N=(3, 3) $nmGWS
	$nmGWS[0][0] = "Low-Beta (Short)";			$nmGWS[0][1] = "14";				$nmGWS[0][2] = "0.00089,0.9,0.008,2.02,1.06,0,0,0,0";
	$nmGWS[1][0] = "High-Beta (Long)";			$nmGWS[1][1] = "16.6";				$nmGWS[1][2] = "0.00089,0.9,0.008,20.85,3.4,0,0,0,0"
	$nmGWS[2][0] = "Dispersive (very short)";		$nmGWS[2][1] = "5.9";				$nmGWS[2][2] = "0.00089,0.9,0.008,25.0603,15.9125,6.25476,-0.579367,0.4205,-0.105"

	redimension/N=(6,5) $nmGWBM
	$nmGWBM[0][0] = "No";											$nmGWBM[0][1] = "";							$nmGWBM[0][2] = "";						$nmGWBM[0][3] = "";		$nmGWBM[0][4] = "";	
	$nmGWBM[1][0] = "Yes: 35 mm Dipoles";							$nmGWBM[1][1] = "BDipole35mm.dat";			$nmGWBM[1][2] = "BDipole35mm.dat";		$nmGWBM[1][3] = "5";	$nmGWBM[1][4] = "1.3";	
	$nmGWBM[2][0] = "Yes: 35 mm Dipole Upstr., 90 mm Downstr.";		$nmGWBM[2][1] = "BDipole35mm.dat";			$nmGWBM[2][2] = "BDipole90mm.dat";		$nmGWBM[2][3] = "5";	$nmGWBM[2][4] = "1.3";	
	$nmGWBM[3][0] = "Yes: 90 mm Dipole Upstr., 35 mm Downstr.";		$nmGWBM[3][1] = "BDipole90mm.dat";			$nmGWBM[3][2] = "BDipole35mm.dat";		$nmGWBM[3][3] = "5";	$nmGWBM[3][4] = "1.3";	
	$nmGWBM[4][0] = "Yes: 90 mm Dipoles";							$nmGWBM[4][1] = "BDipole90mm.dat";			$nmGWBM[4][2] = "BDipole90mm.dat";		$nmGWBM[4][3] = "5";	$nmGWBM[4][4] = "1.3";	
	$nmGWBM[5][0] = "Yes: 35 mm Test Dipoles";						$nmGWBM[5][1] = "BDipole35mmNoseTest.dat";$nmGWBM[5][2] = "BDipole35mmNoseTest.dat";$nmGWBM[5][3] = "1";$nmGWBM[5][4] = "1.3";	

	redimension/N=2 $nmGWMP
	$nmGWMP[0] = "NSLSII";						$nmGWMP[1] = "3,0.5,0,0,0,0,0";		//$nmGWMP[2] = ""
endif

SrwUtiSrcElecAndID()
end

//+++++++++++++++++++++++++++++++++++++++
//List of available SR Sources
//+++++++++++++++++++++++++++++++++++++++
function/S SrwUtiSrcListAvail(s)
string s
string strSrcList = "SOLEIL" + s + "NSLS-II" //to edit
return strSrcList
end

//+++++++++++++++++++++++++++++++++++++++
//IDs linked to straight sections (to be called only after SrwUtiSrcSelect(), which defines an SR source!)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiSrcElecAndID(ElecName, MagName, StrSect, WithBM, ID1, s1, ID2, s2, newZeroS, disp)
string ElecName = srwUtiGetValS("SrwElecName", "Elec", "")
string MagName = srwUtiGetValS("SrwMagName", "Mag", "")
variable StrSect = srwUtiGetValN("StrSect", 1, "SrwUtiSrcElecAndID")
variable WithBM = srwUtiGetValN("WithBM", 1, "SrwUtiSrcElecAndID")
variable ID1 =  srwUtiGetValN("ID1", 1, "SrwUtiSrcElecAndID")
variable s1 = srwUtiGetValN("s1", 0., "SrwUtiSrcElecAndID")
variable ID2 =  srwUtiGetValN("ID2", 1, "SrwUtiSrcElecAndID")
variable s2 = srwUtiGetValN("s2", 0., "SrwUtiSrcElecAndID")
variable newZeroS = srwUtiGetValN("newZeroS", 0., "SrwUtiSrcElecAndID")
variable disp = srwUtiGetValN("disp", 1, "SrwUtiSrcElecAndID")
prompt ElecName, "Name for Electron Beam to create"
prompt MagName, "Name for Magnetic Field to create"
prompt StrSect, "Straight Section", popup srwUtiSrcListID("gwSrwUtiSrcStrSect", ";")
prompt WithBM, "Include Bending Magnet Edges", popup srwUtiSrcListID("gwSrwUtiSrcBM", ";")
prompt ID1, "First ID", popup " ;" + srwUtiSrcListID("gwSrwUtiSrcIDs", ";")
prompt s1, "Center Position of the First ID"
prompt ID2, "Second ID", popup " ;" + srwUtiSrcListID("gwSrwUtiSrcIDs", ";")
prompt s2, "Center Position of the Second ID"
prompt newZeroS, "Origin of Longitudinal Coord. [m]"
prompt disp, "Display?", popup "No;Yes: Magnetic Field;Yes: Electron Trajectory;Yes: Magn. Field and Trajectory"
Silent 1						|	Setting up Magnetic Field  ...
PauseUpdate
srwUtiSetValS("SrwElecName", ElecName, "")
srwUtiSetValS("SrwMagName", MagName, "")
srwUtiSetValS("SrwMagGenTotName", MagName + SrwFieldType, "")
srwUtiSetValN("StrSect", StrSect, "SrwUtiSrcElecAndID")
srwUtiSetValN("WithBM", WithBM, "SrwUtiSrcElecAndID")
srwUtiSetValN("ID1", ID1, "SrwUtiSrcElecAndID")
srwUtiSetValN("s1", s1, "SrwUtiSrcElecAndID")
srwUtiSetValN("ID2", ID2, "SrwUtiSrcElecAndID")
srwUtiSetValN("s2", s2, "SrwUtiSrcElecAndID")
srwUtiSetValN("newZeroS", newZeroS, "SrwUtiSrcElecAndID")
srwUtiSetValN("disp", disp, "SrwUtiSrcElecAndID")

variable sNp = 40001 //total number of points in resulting field

string nmGWID = "gwSrwUtiSrcIDs", nmGWS = "gwSrwUtiSrcStrSect", nmGWBM = "gwSrwUtiSrcBM", nmGWMP = "gwSrwUtiSrcMainParam"
if((exists(nmGWID) != 1) %| (exists(nmGWS) != 1) %| (exists(nmGWBM) != 1) %| (exists(nmGWMP) != 1))
	abort "Can not find Source information"
endif

//PathInfo SrwUtiPath
//string AbsPathToMagFieldDir = S_path + "SRW Procedures:Utilities:MagnFieldIDs:" + $nmGWMP[0] + ":"
string AbsPathToMagFieldDir = srwUtiPathStr() + "SRW Procedures:Utilities:MagnFieldIDs:" + $nmGWMP[0] + ":"

make/O/T/N=6 wAuxFieldLoadedUtiSrc
make/O wNumIDUtiSrc = {ID1, ID2}, wPosIDUtiSrc = {s1, s2}

variable sStart = 1E+23, sEnd = -1E+23, sStartCur = 0, sEndCur = 0, sRange = 0, sNpB = 0, sStep = 0
string nmFileB, nmFieldWaveCore = "wAuxFieldSrwUtiSrc", nmFieldWave
variable i = 0, curID, curPos, iComp
do
	curID = wNumIDUtiSrc[i]
	curPos = wPosIDUtiSrc[i]
	
	if(curID > 1)
		sStep = str2num($nmGWID[curID - 2][3])
	
		iComp = 1
		do
			nmFileB = $nmGWID[curID - 2][iComp]
			if(strlen(nmFileB) > 0)
				nmFileB = AbsPathToMagFieldDir + nmFileB
				LoadWave/G/D/O/A=wMagImportAux/N/Q nmFileB
				if(V_flag != 1)
					abort "Failed reading magnetic field data file"
				endif
				nmFieldWave = nmFieldWaveCore + num2str(i) + num2str(iComp)
				duplicate/O wMagImportAux0 $nmFieldWave
				killwaves/Z wMagImportAux0
				 
				wAuxFieldLoadedUtiSrc[i*2 + iComp - 1] = nmFieldWave
				sNpB = dimsize($nmFieldWave, 0)
				sRange = (sNpB - 1)*sStep*0.001
				sStartCur = curPos - 0.5*sRange
				sEndCur = curPos + 0.5*sRange
				SetScale/I x sStartCur, sEndCur, "m", $nmFieldWave
				
				if(sStart > sStartCur)
					sStart = sStartCur
				endif
				if(sEnd < sEndCur)
					sEnd = sEndCur
				endif
			endif
			iComp += 1
		while(iComp <= 2)
	endif
	i += 1
while(i < 2)

//Setting-up electron beam
string elecNameTot = ElecName + "_ebm"
string str2exe = "SrwElecFilament(\"" + ElecName + "\"," + $nmGWMP[1] + ");"
str2exe += "SrwElecThick(\"" + elecNameTot + "\"," + $nmGWS[StrSect - 1][2] + ")"
execute/Q str2exe

variable sStep_m = sStep*0.001
variable s0ElecBeam = sStart + sStep_m
string nmCoreAuxDrift = "AuxDriftUtiSrc"
string nmAuxDrift = nmCoreAuxDrift + "_bli"
if((sStart < 1E+23) %& (abs(s0ElecBeam) > sStep_m))
	SrwOptDrift(nmCoreAuxDrift, s0ElecBeam) //assuming default position s0 = 0
	srElecBeamPropag($elecNameTot, $nmAuxDrift)
endif

variable strSectLen = str2num($nmGWS[StrSect - 1][1])
variable sCen, sAux
variable halfDipoleIronLen = 0
if(WithBM > 1)
	sStep = str2num($nmGWBM[WithBM - 1][3])
	halfDipoleIronLen = str2num($nmGWBM[WithBM - 1][4])

	nmFileB = $nmGWBM[WithBM - 1][1] //upstream
	if(strlen(nmFileB) > 0)
		nmFileB = AbsPathToMagFieldDir + nmFileB
		LoadWave/G/D/O/A=wMagImportAux/N/Q nmFileB
		if(V_flag != 1)
			abort "Failed reading magnetic field data file"
		endif
		nmFieldWave = nmFieldWaveCore + "DU"
		duplicate/O wMagImportAux0 $nmFieldWave
		killwaves/Z wMagImportAux0

		wAuxFieldLoadedUtiSrc[4] = nmFieldWave
		sNpB = dimsize($nmFieldWave, 0)
		sRange = (sNpB - 1)*sStep*0.001
		sCen = -0.5*strSectLen - halfDipoleIronLen
		sStartCur = sCen - 0.5*sRange
		sEndCur = sCen + 0.5*sRange
		SetScale/I x sStartCur, sEndCur, "m", $nmFieldWave
		sAux = sCen + 0.5*halfDipoleIronLen
		if(sStart > sAux)
			sStart = sAux
		endif
	endif

	nmFileB = $nmGWBM[WithBM - 1][2] //upstream
	if(strlen(nmFileB) > 0)
		nmFileB = AbsPathToMagFieldDir + nmFileB
		LoadWave/G/D/O/A=wMagImportAux/N/Q nmFileB
		if(V_flag != 1)
			abort "Failed reading magnetic field data file"
		endif
		nmFieldWave = nmFieldWaveCore + "DD"
		duplicate/O wMagImportAux0 $nmFieldWave
		killwaves/Z wMagImportAux0

		wAuxFieldLoadedUtiSrc[5] = nmFieldWave
		sNpB = dimsize($nmFieldWave, 0)
		sRange = (sNpB - 1)*sStep*0.001
		sCen = 0.5*strSectLen + halfDipoleIronLen
		sStartCur = sCen - 0.5*sRange
		sEndCur = sCen + 0.5*sRange
		SetScale/I x sStartCur, sEndCur, "m", $nmFieldWave
		sAux = sCen - 0.5*halfDipoleIronLen
		if(sEnd < sAux)
			sEnd = sAux
		endif
	endif
endif

SrwMagFieldCreate(MagName, 0.5*(sStart + sEnd), sEnd - sStart, sNp)
string nmMagBX = MagName + "BX_fld", nmMagBZ = MagName + "BZ_fld", magNameTot = MagName + "_mag"
$nmMagBX = 0; $nmMagBZ = 0
i = 0
do
	nmFieldWave = wAuxFieldLoadedUtiSrc[i]
	if(exists(nmFieldWave) == 1)
		sStartCur = dimoffset($nmFieldWave, 0); sEndCur = sStartCur + (dimsize($nmFieldWave, 0) - 1)*dimdelta($nmFieldWave, 0)
		if((i == 0) %| (i == 2))
			$nmMagBX += $nmFieldWave(x)*srwUtiNonZeroIntervB(x, sStartCur, sEndCur)
		else
			$nmMagBZ += $nmFieldWave(x)*srwUtiNonZeroIntervB(x, sStartCur, sEndCur)
		endif
		killwaves/Z $nmFieldWave
	endif
	i += 1
while(i < 6)

variable sStepFin = dimdelta($nmMagBZ, 0)
variable sStartOld
if(abs(newZeroS) > abs(sStepFin))
	sStartOld = dimoffset($nmMagBZ, 0)
	SetScale/P x (sStartOld - newZeroS), sStepFin, "m", $nmMagBX, $nmMagBZ
	s0ElecBeam = srwGetElecBeamLongPos(elecNameTot)
	srwSetElecBeamLongPos(elecNameTot, s0ElecBeam - newZeroS)
endif

if((disp == 2) %| (disp == 4)) //display magnetic field
	SrwMagDisplayField(nmMagBX); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(0,10,350,200,0,0)
	SrwMagDisplayField(nmMagBZ); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(200,10,350,200,0,0)
endif
string nmCoreTraj = srwUtiTruncString(ElecName + MagName, 26)
if((disp == 3) %| (disp == 4)) //display trajectory
	SrwTrjCreateTransvUnif(nmCoreTraj, elecNameTot, magNameTot, 2,1,1,1); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(0,210,350,200,0,0)
	SrwTrjCreateTransvUnif(nmCoreTraj, elecNameTot, magNameTot, 1,2,1,1); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(150,210,350,200,0,0)
	SrwTrjCreateTransvUnif(nmCoreTraj, elecNameTot, magNameTot, 1,1,2,1); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(300,210,350,200,0,0)
	SrwTrjCreateTransvUnif(nmCoreTraj, elecNameTot, magNameTot, 1,1,1,2); SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,210,350,200,0,0)
endif

killwaves/Z wAuxFieldLoadedUtiSrc, wNumIDUtiSrc, wPosIDUtiSrc //, $nmGWID, $nmGWS, $nmGWBM, $nmGWMP
end

//+++++++++++++++++++++++++++++++++++++++
//List of Straight Sections or Insertion Devices available at some Sources
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiSrcListID(nmWave, s)
string nmWave, s
string strRes = ""

//string nmGWID = "gwSrwUtiSrcIDs"
if(exists(nmWave) != 1)
	return strRes
endif
wave/T wID = $nmWave

variable nID = dimsize(wID, 0)
if(nID <= 0)
	return strRes
endif
variable i, nIDm1 = nID - 1
for(i = 0; i < nID; i += 1)
	strRes += wID[i][0]
	if(i < nIDm1)
		strRes += s
	endif
endfor
return strRes
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete: Soleil IDs linked to straight sections
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiSrcElecAndID_SOLEIL(ElecName, MagName, StrSect, WithBM, ID1, s1, ID2, s2, ID3, s3)
string ElecName = srwUtiGetValS("SrwElecName", "Elec", "")
string MagName = srwUtiGetValS("SrwMagName", "Mag", "")
variable StrSect = srwUtiGetValN("StrSect", 2, "SolElecAndMagInStrSect")
variable WithBM = srwUtiGetValN("WithBM", 1, "SolElecAndMagInStrSect")
string ID1 =  srwUtiGetValS("ID1", "", "SolElecAndMagInStrSect")
variable s1 = srwUtiGetValN("s1", 1.85, "SolElecAndMagInStrSect")
string ID2 =  srwUtiGetValS("ID2", "", "SolElecAndMagInStrSect")
variable s2 = srwUtiGetValN("s2", 0., "SolElecAndMagInStrSect")
string ID3 =  srwUtiGetValS("ID3", "", "SolElecAndMagInStrSect")
variable s3 = srwUtiGetValN("s3", -1.85, "SolElecAndMagInStrSect")
prompt ElecName, "Name of Electron Beam to create"
prompt MagName, "Name of Magnetic Field to create"
prompt StrSect, "Straight Section", popup "Short;Medium;Long;Very Short"
prompt WithBM, "Include Bending Magnet Edges", popup"No;Yes"
prompt ID1, "First ID", popup " ;" + srwUtiSrcIDNamesList_SOLEIL(";")
prompt s1, "Center Position of the First ID"
prompt ID2, "Second ID", popup " ;" + srwUtiSrcIDNamesList_SOLEIL(";")
prompt s2, "Center Position of the Second ID"
prompt ID3, "Third ID", popup " ;" + srwUtiSrcIDNamesList_SOLEIL(";")
prompt s3, "Center Position of the Third ID"
Silent 1						|	Setting up Magnetic Field  ...
PauseUpdate
srwUtiSetValS("SrwElecName", ElecName, "")
srwUtiSetValS("SrwMagName", MagName, "")
srwUtiSetValS("SrwMagGenTotName", MagName + SrwFieldType, "")
srwUtiSetValN("StrSect", StrSect, "SolElecAndMagInStrSect")
srwUtiSetValN("WithBM", WithBM, "SolElecAndMagInStrSect")
srwUtiSetValS("ID1", ID1, "SolElecAndMagInStrSect")
srwUtiSetValN("s1", s1, "SolElecAndMagInStrSect")
srwUtiSetValS("ID2", ID2, "SolElecAndMagInStrSect")
srwUtiSetValN("s2", s2, "SolElecAndMagInStrSect")
srwUtiSetValS("ID3", ID3, "SolElecAndMagInStrSect")
srwUtiSetValN("s3", s3, "SolElecAndMagInStrSect")

variable VeryShortSectLen = 4.4527 //[m] //no room for ID in this section

variable ShortSectLen = 7.722 //[m]
variable ShortSectLenBwQuad = 3.649 //[m]

variable MedSectLen = 12.405 //[m]
variable MedSectLenBwQuad = 6.964 //[m]

variable LongSectLen = 18.54 //[m]
variable LongSectLenBwQuad = 11.8 //[m]

variable BMEdgeLen = 80. //[mm]
variable BMConstField =1.718 //[T]
variable BMConstFieldLen = 0.9 //0.7//0.4 //[m]
string PathToMagFiledFiles = "SRW Procedures:Utilities:MagnFieldIDs:SOLEIL:"

variable MagFieldRange = 0, DistBwBMEdges = 0, DistBwQuads = 0

SrwUtiTriggerPrint(2)
SrwElecFilament(ElecName,2.75,0.5,0,0,0,0,0)
string TotElecName = ElecName + "_ebm"

if(StrSect == 1) //Short
	SrwElecThick(TotElecName,0.001016,3.73,0.037,17.78,1.75,0,0,0.28,0)
	
	DistBwQuads = ShortSectLenBwQuad
	MagFieldRange = ShortSectLenBwQuad
	if(WithBM == 2)
		DistBwBMEdges = ShortSectLen //- 0.002*BMEdgeLen
		MagFieldRange = ShortSectLen + 2*BMConstFieldLen
	endif
endif
if(StrSect == 2) //Medium
	SrwElecThick(TotElecName,0.001016,3.73,0.037,4,1.77,0,0,0.13,0)

	DistBwQuads = MedSectLenBwQuad
	MagFieldRange = MedSectLenBwQuad
	if(WithBM == 2)
		DistBwBMEdges = MedSectLen //- 0.002*BMEdgeLen
		MagFieldRange = MedSectLen + 2*BMConstFieldLen
	endif
endif
if(StrSect == 3) //Long
	SrwElecThick(TotElecName,0.001016,3.73,0.037,10.09,8.01,0,0,0.2,0)

	DistBwQuads = LongSectLenBwQuad
	MagFieldRange = LongSectLenBwQuad
	if(WithBM == 2)
		DistBwBMEdges = LongSectLen //- 0.002*BMEdgeLen
		MagFieldRange = LongSectLen + 2*BMConstFieldLen
	endif
endif
if(StrSect == 4) //Very Short
	SrwElecThick(TotElecName,0.001016,3.73,0.037,17.78,1.75,0,0,0.28,0)
	//the data taken from Short section

	DistBwQuads = 0
	MagFieldRange = 0
	if(WithBM == 2)
		DistBwBMEdges = VeryShortSectLen
		MagFieldRange = VeryShortSectLen + 2*BMConstFieldLen
	endif
endif

//PathInfo igor
PathInfo SrwUtiPath
string AbsPathToMagFieldDir = S_path + PathToMagFiledFiles
SolIDCreateNames()

string CurSegmName, CurFileNameX, CurFileNameZ, CurAbsPath

SrwMagFieldCreate(MagName,0,MagFieldRange,20000)
SrwMagZero(MagName + "BX_fld"); SrwMagZero(MagName + "BZ_fld")
if(WithBM == 2)
	//SrwMagEdge(MagName + "BZ_fld",2,0,DistBwBMEdges,BMEdgeLen,1.718)
	
	variable DipoleMagDataStep = 5 //1 //[mm]
	variable DipoleHalfIronLength = 0.52537 //0.5180 //[m]
	CurFileNameZ = "BzSolDipole.dta"
	
	SrwMagAddImpCmpn(MagName + "BZ_fld",2,0.5*DistBwBMEdges + DipoleHalfIronLength,DipoleMagDataStep,1,AbsPathToMagFieldDir + CurFileNameZ)
	SrwMagAddImpCmpn(MagName + "BZ_fld",2,-0.5*DistBwBMEdges - DipoleHalfIronLength,DipoleMagDataStep,1,AbsPathToMagFieldDir + CurFileNameZ)
endif

if(cmpstr(ID1, " ") == 0)
	ID1 = ""
endif
if(cmpstr(ID2, " ") == 0)
	ID2 = ""
endif
if(cmpstr(ID3, " ") == 0)
	ID3 = ""
endif

variable AmOfSegm = 0
make/O/T/N=3 AuxID
make/O/N=3 AuxCenPos
if(strlen(ID1) > 0)
	AuxID[AmOfSegm] = ID1
	AuxCenPos[AmOfSegm] = s1
	AmOfSegm += 1
endif
if(strlen(ID2) > 0)
	AuxID[AmOfSegm] = ID2
	AuxCenPos[AmOfSegm] = s2
	AmOfSegm += 1
endif
if(strlen(ID3) > 0)
	AuxID[AmOfSegm] = ID3
	AuxCenPos[AmOfSegm] = s3
	AmOfSegm += 1
endif

if(AmOfSegm > 0)
	string AuxDriftName = "Drift"
	string TotAuxDriftName = AuxDriftName + "_bli"
	SrwOptDrift("Drift",-0.499*DistBwQuads)
	srElecBeamPropag($TotElecName, $TotAuxDriftName)
endif

variable nSegm = 0
variable CurSegmPos, CurStep
string MagFieldDataPath
do
	CurSegmName = AuxID[nSegm]
	CurSegmPos = AuxCenPos[nSegm]

	if(cmpstr(CurSegmName, SolNameU20) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_U20.dta"
		CurStep = 0.2 //0.5
	endif
	if(cmpstr(CurSegmName, SolNameHU40LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU40_LinHor.dta"
		CurStep = 0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU40HE) == 0)
		CurFileNameX = "Bx_HU40_Hel.dta"
		CurFileNameZ = "Bz_HU40_Hel.dta"
		CurStep = 0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU40LV) == 0)
		CurFileNameX = "Bx_HU40_LinVert.dta"
		CurFileNameZ = "Bz_HU40_LinVert.dta"
		CurStep = 0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU45LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU45_LinHor.dta"
		CurStep = 0.9225
	endif
	if(cmpstr(CurSegmName, SolNameHU45HE) == 0)
		CurFileNameX = "Bx_HU45_Hel.dta"
		CurFileNameZ = "Bz_HU45_Hel.dta"
		CurStep = 0.9225
	endif
	if(cmpstr(CurSegmName, SolNameHU45LV) == 0)
		CurFileNameX = "Bx_HU45_LinVert.dta"
		CurFileNameZ = "Bz_HU45_LinVert.dta"
		CurStep = 0.9225
	endif
	if(cmpstr(CurSegmName, SolNameHU60LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU60_15_LinHor.dta"
		CurStep = 1.05//0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU60HE) == 0)
		CurFileNameX = "Bx_HU60_15_Hel.dta"
		CurFileNameZ = "Bz_HU60_15_Hel.dta"
		CurStep = 1.05//0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU60LV) == 0)
		CurFileNameX = "Bx_HU60_15_LinVert.dta"
		CurFileNameZ = "Bz_HU60_15_LinVert.dta"
		CurStep = 1.05//0.94
	endif
	if(cmpstr(CurSegmName, SolNameHU70LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU70_LinHor.dta"
		CurStep = 1.015
	endif
	if(cmpstr(CurSegmName, SolNameHU70HE) == 0)
		CurFileNameX = "Bx_HU70_Hel.dta"
		CurFileNameZ = "Bz_HU70_Hel.dta"
		CurStep = 1.015
	endif
	if(cmpstr(CurSegmName, SolNameHU70LV) == 0)
		CurFileNameX = "Bx_HU70_LinVert.dta"
		CurFileNameZ = "Bz_HU70_LinVert.dta"
		CurStep = 1.015
	endif
	if(cmpstr(CurSegmName, SolNameHU80LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU80_ELETTRA_15_LinHor.dta" //"Bz_HU80_LinHor.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameHU80HE) == 0)
		CurFileNameX = "Bx_HU80_ELETTRA_15_Hel.dta" //"Bx_HU80_Hel.dta"
		CurFileNameZ = "Bz_HU80_ELETTRA_15_Hel.dta" //"Bz_HU80_Hel.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameHU80LV) == 0)
		CurFileNameX = "Bx_HU80_ELETTRA_15_LinVert.dta" //"Bx_HU80_LinVert.dta"
		CurFileNameZ = "Bz_HU80_ELETTRA_15_LinVert.dta" //"Bz_HU80_LinVert.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameHU80LH18) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU80_18_LinHor.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameHU80HE18) == 0)
		CurFileNameX = "Bx_HU80_18_Hel.dta"
		CurFileNameZ = "Bz_HU80_18_Hel.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameHU80LV18) == 0)
		CurFileNameX = "Bx_HU80_18_LinVert.dta"
		CurFileNameZ = "Bz_HU80_18_LinVert.dta"
		CurStep = 1.08
	endif
	if(cmpstr(CurSegmName, SolNameW100) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_W100.dta"
		CurStep = 1.
	endif
	if(cmpstr(CurSegmName, SolNameHU256LH) == 0)
		CurFileNameX = "Bx_HU256_0A2650A.dta"
		CurFileNameZ = "Bz_HU256_0A2650A.dta"
		CurStep = 1.28
	endif
	if(cmpstr(CurSegmName, SolNameHU256HE) == 0)
		CurFileNameX = "Bx_HU256_Hel.dta"  //"Bx_HU256_2718A845A.dta"
		CurFileNameZ = "Bz_HU256_Hel.dta"  //"Bz_HU256_2718A845A.dta"
		CurStep = 1.28
	endif
	if(cmpstr(CurSegmName, SolNameHU256LV) == 0)
		CurFileNameX = "Bx_HU256_Hel.dta" //"Bx_HU256_6000A0A.dta"
		CurFileNameZ = ""
		CurStep = 1.28
	endif
	if(cmpstr(CurSegmName, SolNameHU600LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU600_Plan.dta"
		CurStep = 5.
	endif
	if(cmpstr(CurSegmName, SolNameHU600HE) == 0)
		CurFileNameX = "Bx_HU600_Hel.dta"
		CurFileNameZ = "Bz_HU600_Hel.dta"
		CurStep = 5.
	endif
	if(cmpstr(CurSegmName, SolNameHU600LV) == 0)
		CurFileNameX = "Bx_HU600_Plan.dta"
		CurFileNameZ = ""
		CurStep = 5.
	endif
	if(cmpstr(CurSegmName, SolNameHU600LT) == 0)
		CurFileNameX = "Bx_HU600_Tilt.dta"
		CurFileNameZ = "Bz_HU600_Tilt.dta"
		CurStep = 5.
	endif
	if(cmpstr(CurSegmName, SolNameHU640LH) == 0)
		CurFileNameX = ""
		CurFileNameZ = "Bz_HU640_Plan.dta"
		CurStep = 6.4
	endif
	if(cmpstr(CurSegmName, SolNameHU640HE) == 0)
		CurFileNameX = "Bx_HU640_Plan.dta"
		CurFileNameZ = "Bz_HU640_Plan.dta"
		CurStep = 6.4
	endif
	if(cmpstr(CurSegmName, SolNameHU640LV) == 0)
		CurFileNameX = "Bx_HU640_Plan.dta"
		CurFileNameZ = ""
		CurStep = 6.4
	endif
	if(cmpstr(CurSegmName, SolNameHU640LT) == 0)
		CurFileNameX = "Bx_HU640_Plan.dta"
		CurFileNameZ = "Bz_HU640_Tilt.dta"
		CurStep = 6.4
	endif
	
	if(strlen(CurFileNameX) > 0)
		SrwMagAddImpCmpn(MagName + "BX_fld",2,CurSegmPos,CurStep,1,AbsPathToMagFieldDir + CurFileNameX)
	endif
	if(strlen(CurFileNameZ) > 0)
		SrwMagAddImpCmpn(MagName + "BZ_fld",2,CurSegmPos,CurStep,1,AbsPathToMagFieldDir + CurFileNameZ)
	endif
	
	nSegm += 1
while(nSegm < AmOfSegm)
KillWaves/Z AuxID, AuxCenPos

SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//(Re-)Creates names of Soleil IDs
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiSrcIDNames_SOLEIL()
string/G SolNameU20 = "U20"
string/G SolNameHU45LH = "HU45 Linear Horizontal"
string/G SolNameHU45HE = "HU45 Helical"
string/G SolNameHU45LV = "HU45 Linear Vertical"
string/G SolNameHU40LH = "HU40 Linear Horizontal"
string/G SolNameHU40HE = "HU40 Helical"
string/G SolNameHU40LV = "HU40 Linear Vertical"
string/G SolNameHU60LH = "HU60 Linear Horizontal"
string/G SolNameHU60HE = "HU60 Helical"
string/G SolNameHU60LV = "HU60 Linear Vertical"
string/G SolNameHU70LH = "HU70 Linear Horizontal"
string/G SolNameHU70HE = "HU70 Helical"
string/G SolNameHU70LV = "HU70 Linear Vertical"
string/G SolNameHU80LH = "HU80 Linear Horizontal"
string/G SolNameHU80HE = "HU80 Helical"
string/G SolNameHU80LV = "HU80 Linear Vertical"
string/G SolNameHU80LH18 = "HU80 Gap 18 Lin. Hor."
string/G SolNameHU80HE18 = "HU80 Gap 18 Helical"
string/G SolNameHU80LV18 = "HU80 Gap 18 Lin. Vert."
string/G SolNameW100 = "W100"
string/G SolNameHU256LH = "HU256 Linear Horizontal"
string/G SolNameHU256HE = "HU256 Helical"
string/G SolNameHU256LV = "HU256 Linear Vertical"
string/G SolNameHU600LH = "HU600 Linear Horizontal"
string/G SolNameHU600HE = "HU600 Helical"
string/G SolNameHU600LV = "HU600 Linear Vertical"
string/G SolNameHU600LT = "HU600 Linear Tilted"
string/G SolNameHU640LH = "HU640 Linear Horizontal"
string/G SolNameHU640HE = "HU640 Helical"
string/G SolNameHU640LV = "HU640 Linear Vertical"
string/G SolNameHU640LT = "HU640 Linear Tilted"

end

//+++++++++++++++++++++++++++++++++++++++
//Lists code names of all Soleil IDs
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiSrcIDNamesList_SOLEIL(s)
string s
execute "SrwUtiSrcIDNames_SOLEIL()"
SVAR SolNameU20, SolNameHU45LH, SolNameHU45HE, SolNameHU45LV, SolNameHU70LH, SolNameHU70HE, SolNameHU70LV
SVAR SolNameHU80LH, SolNameHU80HE, SolNameHU80LV, SolNameHU80LH18, SolNameHU80HE18, SolNameHU80LV18
SVAR SolNameW100, SolNameHU256LH, SolNameHU256HE, SolNameHU256LV
SVAR SolNameHU600LH, SolNameHU600HE, SolNameHU600LV, SolNameHU600LT
SVAR SolNameHU640LH, SolNameHU640HE, SolNameHU640LV, SolNameHU640LT
SVAR SolNameHU40LH, SolNameHU40HE, SolNameHU40LV, SolNameHU60LH, SolNameHU60HE, SolNameHU60LV

string AllNames = SolNameU20 + s 
AllNames += SolNameHU40LH + s + SolNameHU40HE + s + SolNameHU40LV + s
AllNames += SolNameHU45LH + s + SolNameHU45HE + s + SolNameHU45LV + s
AllNames += SolNameHU60LH + s + SolNameHU60HE + s + SolNameHU60LV + s
AllNames += SolNameHU70LH + s + SolNameHU70HE + s + SolNameHU70LV + s
AllNames += SolNameHU80LH + s + SolNameHU80HE + s + SolNameHU80LV + s
AllNames += SolNameHU80LH18 + s + SolNameHU80HE18 + s + SolNameHU80LV18 + s
AllNames += SolNameW100 + s + SolNameHU256LH + s + SolNameHU256HE + s + SolNameHU256LV + s
AllNames += SolNameHU600LH + s + SolNameHU600HE + s + SolNameHU600LV + s + SolNameHU600LT + s
AllNames += SolNameHU640LH + s + SolNameHU640HE + s + SolNameHU640LV + s + SolNameHU640LT

return AllNames
end
