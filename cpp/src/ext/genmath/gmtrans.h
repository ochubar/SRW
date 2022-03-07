/************************************************************************//**
 * File: gmtrans.h
 * Description: Space transformations / symmetries (header)
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMTRANS_H
#define __GMTRANS_H

#include "gmvect.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class gmTrans {
protected:
	TMatrix3d M, M_inv;
	TVector3d V;
	double detM, s;

public:
	int ID_No;

	gmTrans(const TMatrix3d& InM, const TMatrix3d& InM_inv,
			const TVector3d& InV, double In_detM, double In_s, int InID_No =-1)
	{
		M = InM; M_inv = InM_inv; V = InV; s = In_s; detM = In_detM; ID_No = InID_No;
	}
	gmTrans(const TMatrix3d& InM, const TVector3d& InV, double In_detM, double In_s, int InID_No =-1)
	{
		M = InM; V = InV; s = In_s; detM = In_detM; M_inv = Matrix3d_inv(InM); ID_No = InID_No;
	}
	gmTrans(TMatrix3d& InM, const TVector3d& InV, double In_s =1, int InID_No =-1)
	{
		M = InM; V = InV; s = In_s; detM = InM.det(); M_inv = Matrix3d_inv(InM); ID_No = InID_No;
	}
	gmTrans() {}
	virtual ~gmTrans() {} //GCC 4.2 produces warning without this

	//int Type_g() { return 2;}
	virtual int Type_Trans() { return 0;}

	//void Dump(std::ostream& o, int ShortSign) // Porting
	//{
	//	radTg::Dump(o);
	//	o << "Transformation: ";
	//	if(ID_No == 1) o << "Translation";
	//	else if(ID_No == 2) o << "Rotation";
	//	else if(ID_No == 3) o << "Plane symmetry";
	//	else if(ID_No == 4) o << "Field inversion";
	//	else if(ID_No == 10) o << "Composite";
	//	if(ShortSign) return;
	//	o << endl;
	//	o << "   Memory occupied: " << SizeOfThis() << " bytes";
	//}

	virtual TVector3d TrPoint(const TVector3d& P) { return M*P+V;}
	virtual TVector3d TrBiPoint(const TVector3d& P) { return M*P;}
	virtual TVector3d TrVectField(const TVector3d& B) { return s*(M*B);}
	virtual TVector3d TrVectPoten(const TVector3d& A) { return /* s*detM*(M_inv*A); */ s*detM*(M*A);}
	virtual TVector3d TrPoint_inv(const TVector3d& P) { return M_inv*(P-V);}
	virtual TVector3d TrBiPoint_inv(const TVector3d& P) { return M_inv*P;}
	virtual TVector3d TrVectField_inv(const TVector3d& B) { return s*(M_inv*B);}
	virtual TVector3d TrVectPoten_inv(const TVector3d& A) { return /* (s/detM)*(M*A); */ (s/detM)*(M_inv*A);}

	virtual TVector3d TrAxialVect(const TVector3d& A) { return detM*(M*A);}
	virtual TVector3d TrAxialVect_inv(const TVector3d& A) { return (1./detM)*(M_inv*A);}

	virtual void TrMatrixGeom(TMatrix3d& Matrix) { Matrix = M*Matrix;}
	virtual void TrMatrixGeom_inv(TMatrix3d& Matrix) { Matrix = M_inv*Matrix;}
	virtual void TrMatrix(TMatrix3d& Matrix) { Matrix = s*M*Matrix;}
	virtual void TrMatrix_inv(TMatrix3d& Matrix) { Matrix = s*M_inv*Matrix;}

	virtual void TrMatrixGeomLeft(TMatrix3d& Matrix) { Matrix = Matrix*M;}
	virtual void TrMatrixGeomLeft_inv(TMatrix3d& Matrix) { Matrix = Matrix*M_inv;}
	virtual void TrMatrixLeft(TMatrix3d& Matrix) { Matrix = s*Matrix*M;}
	virtual void TrMatrixLeft_inv(TMatrix3d& Matrix) { Matrix = s*Matrix*M_inv;}

	int ShowParity() { return int(detM);}
	TMatrix3d& OutMatrix() { return M;}
	TVector3d& OutVector() { return V;}

	TMatrix3d& GetMatrix() { return M;} 
	TMatrix3d& GetMatrix_inv() { return M_inv;} 
	TVector3d& GetVector() { return V;} 
	TVector3d GetVector_inv() 
	{ 
		TVector3d Zero(0.,0.,0.);
		return M_inv*(Zero-V);
	} 

	void SetVector(const TVector3d& InV) { V = InV;} 
	void SetMatrixVector(TMatrix3d& InM, const TVector3d& InV, double In_s =1, int InID_No =-1)
	{
		M = InM; V = InV; s = In_s; detM = InM.det(); M_inv = Matrix3d_inv(InM); ID_No = InID_No;
	}

	//radTField TrField(const radTField& InField)
	//{
	//	radTField OutField(InField);
	//	if(InField.FieldKey.B_) OutField.B = TrVectField(InField.B);
	//	if(InField.FieldKey.H_) OutField.H = TrVectField(InField.H);
	//	if(InField.FieldKey.A_) OutField.A = TrVectPoten(InField.A);
	//	if(InField.FieldKey.M_) OutField.M = TrVectField(InField.M);
	//	if(InField.FieldKey.Ib_) OutField.Ib = TrVectField(InField.Ib);
	//	if(InField.FieldKey.Ih_) OutField.Ih = TrVectField(InField.Ih);
	//	return OutField;
	//}
	//radTField TrField_inv(const radTField& InField)
	//{
	//	radTField OutField(InField);
	//	if(InField.FieldKey.B_) OutField.B = TrVectField_inv(InField.B);
	//	if(InField.FieldKey.H_) OutField.H = TrVectField_inv(InField.H);
	//	if(InField.FieldKey.A_) OutField.A = TrVectPoten_inv(InField.A);
	//	if(InField.FieldKey.M_) OutField.M = TrVectField_inv(InField.M);    // Really _inv here ?
	//	if(InField.FieldKey.Ib_) OutField.Ib = TrVectField_inv(InField.Ib);
	//	if(InField.FieldKey.Ih_) OutField.Ih = TrVectField_inv(InField.Ih);
	//	return OutField;
	//}
	void Invert()
	{
		TMatrix3d OldM = M;
		M = M_inv; M_inv = OldM;
		V = (-1.)*(M*V);
		detM = 1./detM;
	}

	//int SizeOfThis() { return sizeof(gmTrans);}

	friend gmTrans Product(const gmTrans& Tr2, const gmTrans& Tr1);
	friend void TrProduct(gmTrans* Tr2Ptr, gmTrans* Tr1Ptr, gmTrans& ResTrPtr);
	
	void TrMult(gmTrans& multTr, char Left_or_Right = 'L')
	{
		if(Left_or_Right == 'L')
		{
			M = multTr.M*M;
			M_inv = M_inv*multTr.M_inv;
			V = multTr.M*V + multTr.V;
			detM *= multTr.detM;
			s *= multTr.s;
		}
		else
		{
			V += M*multTr.V;
			M = M*multTr.M;
			M_inv = multTr.M_inv*M_inv;
			detM *= multTr.detM;
			s *= multTr.s;
		}
		//ResTr.M = Tr2Ptr->M*Tr1Ptr->M;
		//ResTr.M_inv = Tr1Ptr->M_inv*Tr2Ptr->M_inv;
		//ResTr.V = Tr2Ptr->M*Tr1Ptr->V+Tr2Ptr->V;
		//ResTr.detM = Tr2Ptr->detM*Tr1Ptr->detM;
		//ResTr.s = Tr2Ptr->s*Tr1Ptr->s;

		ID_No = 10; // Composite
	}

	void SetupRotation(const TVector3d&, const TVector3d&, double);
	void SetupRotation(const TVector3d& vPointOnAxis, const TVector3d& vUnit1, const TVector3d& vUnit2);
	void SetupRotationToPermutAxes(const TVector3d& vCenPoint, char DefOrient, char Orient);
	void SetupPlaneSym(const TVector3d& PoiOnPlaneVect, const TVector3d& N);

	void SetupTranslation(const TVector3d& InV)
	{
		TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.);
		//M = TMatrix3d(St0, St1, St2); M_inv = M; V = InV; //causes error depending in struct alignment !!!
		M.Str0 = St0; M.Str1 = St1; M.Str2 = St2;

		//M.Str0.x = 1; M.Str0.y = 0; M.Str0.z = 0; 
		//M.Str0.x = 0; M.Str0.y = 1; M.Str0.z = 0; 
		//M.Str0.x = 0; M.Str0.y = 0; M.Str0.z = 1; 

		M_inv = M; 
		V = InV;
		detM = s = 1.;
		ID_No = 1;
	}
	void SetupIdent()
	{
		TVector3d St0(1.,0.,0.), St1(0.,1.,0.), St2(0.,0.,1.), Zero(0.,0.,0.);
		//M = TMatrix3d(St0, St1, St2); M_inv = M; //causes error depending in struct alignment !!!
		M.Str0 = St0; M.Str1 = St1; M.Str2 = St2;
		M_inv = M;
		V = Zero;
		detM = s = 1.;
		ID_No = 10;
	}

	char IsIdent(double RelTol)
	{
		if((fabs(M.Str0.x - 1.) > RelTol) || (fabs(M.Str0.y) > RelTol) || (fabs(M.Str0.z) > RelTol)) return 0;
		if((fabs(M.Str1.y - 1.) > RelTol) || (fabs(M.Str1.x) > RelTol) || (fabs(M.Str1.z) > RelTol)) return 0;
		if((fabs(M.Str2.z - 1.) > RelTol) || (fabs(M.Str2.x) > RelTol) || (fabs(M.Str2.y) > RelTol)) return 0;
		if((fabs(V.x) > RelTol) || (fabs(V.y) > RelTol) || (fabs(V.z) > RelTol)) return 0;
		if(fabs(detM - 1.) > RelTol) return 0;
		return 1;
	}
};

//-------------------------------------------------------------------------

inline gmTrans Product(const gmTrans& Tr2, const gmTrans& Tr1)
{
	return gmTrans(Tr2.M*Tr1.M, Tr1.M_inv*Tr2.M_inv, Tr2.M*Tr1.V+Tr2.V, Tr2.detM*Tr1.detM, Tr2.s*Tr1.s, 10); // ID_No = 10; // Composite
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class gmIdentTrans : public gmTrans {
public:
	gmIdentTrans() 
	{
		M.Str0.x = 1; M.Str0.y = 0; M.Str0.z = 0;
		M.Str1.x = 0; M.Str1.y = 1; M.Str1.z = 0;
		M.Str2.x = 0; M.Str2.y = 0; M.Str2.z = 1;
		M_inv = M; 
		V.x = V.y = V.z = 0; 
		s = 1; 
		detM = 1; 
	}

	TVector3d TrPoint(const TVector3d& P) { return P;}
	TVector3d TrBiPoint(const TVector3d& P) { return P;}
	TVector3d TrVectField(const TVector3d& B) { return B;}
	TVector3d TrVectPoten(const TVector3d& A) { return A;}
	TVector3d TrPoint_inv(const TVector3d& P) { return P;}
	TVector3d TrBiPoint_inv(const TVector3d& P) { return P;}
	TVector3d TrVectField_inv(const TVector3d& B) { return B;}
	TVector3d TrVectPoten_inv(const TVector3d& A) { return A;}
		
	TVector3d TrAxialVect(const TVector3d& A) { return A;}
	TVector3d TrAxialVect_inv(const TVector3d& A) { return A;}

	void TrMatrix(TMatrix3d& Matrix) {}
	void TrMatrix_inv(TMatrix3d& Matrix) {}

	int Type_Trans() { return 1;}
	//int SizeOfThis() { return sizeof(gmIdentTrans);}

	friend void TrProduct(gmTrans* Tr2Ptr, gmTrans* Tr1Ptr, gmTrans& ResTrPtr);
};

//-------------------------------------------------------------------------

inline void TrProduct(gmTrans* Tr2Ptr, gmTrans* Tr1Ptr, gmTrans& ResTr)
{
	gmIdentTrans IdentTr;
	int IdentTrID = IdentTr.Type_Trans();
	int Tr1ID = Tr1Ptr->Type_Trans();
	int Tr2ID = Tr2Ptr->Type_Trans();

	if(Tr2ID==IdentTrID) 
	{
		ResTr = *Tr1Ptr;
	}
	else if(Tr1ID==IdentTrID)
	{
		ResTr = *Tr2Ptr;
	}
	else
	{
		ResTr.M = Tr2Ptr->M*Tr1Ptr->M;
		ResTr.M_inv = Tr1Ptr->M_inv*Tr2Ptr->M_inv;
		ResTr.V = Tr2Ptr->M*Tr1Ptr->V+Tr2Ptr->V;
		ResTr.detM = Tr2Ptr->detM*Tr1Ptr->detM;
		ResTr.s = Tr2Ptr->s*Tr1Ptr->s;
	}

	ResTr.ID_No = 10; // Composite
}

//-------------------------------------------------------------------------

#endif
