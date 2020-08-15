/************************************************************************//**
 * File: sroptzp.cpp
 * Description: Optical element: Zone Plate
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptzp.h"

//*************************************************************************

void srTZonePlate::DefineAttenModulConstants()
{
	m_aModH = m_bModH = m_cModH = m_dModH = 0;
	m_ModH_IsDefined = false;

	if(Nzones <= 1) return;

	if(m_ZoneHeightRatioExtToCen >= 0)
	{
		m_ModH_IsDefined = true;

		double r1 = sqrt(RnMaxe2/Nzones);
		double r2 = sqrt(RnMaxe2*2/Nzones);
		double rNm1 = sqrt(RnMaxe2*(Nzones - 1)/Nzones);
		double x0 = 0.5*(r1 + r2), xe = 0.5*(rNm1 + RnMax);
		double x1 = 0, x2 = 0;
		double h0 = Thickness, he = m_ZoneHeightRatioExtToCen*Thickness;
		double h1 = 0, h2 = 0;

		if((m_ZoneIntermedNum1 > 0) && (m_ZoneHeightRatioIntermedToCen1 > 0))
		{
            double ri1 = sqrt(RnMaxe2*m_ZoneIntermedNum1/Nzones);
            double ri1p1 = sqrt(RnMaxe2*(m_ZoneIntermedNum1 + 1)/Nzones);
			x1 = 0.5*(ri1 + ri1p1);
			if((x1 == x0) || (x1 == xe)) x1 = 0;
			h1 = m_ZoneHeightRatioIntermedToCen1*Thickness;

			if((m_ZoneIntermedNum2 > 0) && (m_ZoneHeightRatioIntermedToCen2 > 0))
			{
				double ri2 = sqrt(RnMaxe2*m_ZoneIntermedNum2/Nzones);
				double ri2p1 = sqrt(RnMaxe2*(m_ZoneIntermedNum2 + 1)/Nzones);
				x2 = 0.5*(ri2 + ri2p1);
				if((x2 == x0) || (x2 == xe) || (x2 == x1)) x2 = 0;
                h2 = m_ZoneHeightRatioIntermedToCen2*Thickness;
			}
		}

		if((x0 != xe) && (x0 != 0) && (xe != 0))
		{
			double x0_mi_xe = (x0 - xe);
			double inv_x0_mi_xe = 1./x0_mi_xe;

			if((x1 != 0) && (x1 != x0) && (x1 != xe))
			{
				double x1_mi_xe = x1 - xe, x0_mi_x1 = x0 - x1;
				double inv_x1_mi_xe = 1./x1_mi_xe, inv_x0_mi_x1 = 1./x0_mi_x1;

				if((x2 != 0) && (x2 != x1) && (x2 != x0) && (x2 != xe))
				{
					double x2_mi_xe = x2 - xe, x0_mi_x2 = x0 - x2, x1_mi_x2 = x1 - x2;
					double CommonMult6 = 1./(x0_mi_x1*x0_mi_x2*x1_mi_x2*x0_mi_xe*x1_mi_xe*x2_mi_xe);

					m_aModH = (-(he*x0_mi_x1*x0_mi_x2*x1_mi_x2) + h2*x0_mi_x1*x0_mi_xe*x1_mi_xe - (h1*x0_mi_x2*x0_mi_xe - h0*x1_mi_x2*x1_mi_xe)*x2_mi_xe)*CommonMult6;  // /((x0_mi_x1)*(x0_mi_x2)*(x1_mi_x2)*(x0_mi_xe)*(x1_mi_xe)*(x2_mi_xe));
					m_bModH = (he*x0_mi_x1*x0_mi_x2*x1_mi_x2*(x0 + x1 + x2) + x2*(h1*x0*x0_mi_x2*(x0 + x2) + h0*x1*(x2*x2 - x1*x1)) - (h1*(x0*x0*x0 - x2*x2*x2) + h0*(x2*x2*x2 - x1*x1*x1))*xe + (h1*x0_mi_x2 - h0*x1_mi_x2)*xe*xe*xe - h2*x0_mi_x1*x0_mi_xe*x1_mi_xe*(x0 + x1 + xe))*CommonMult6; // /((x0 - x1)*(x0 - x2)*(x1 - x2)*(x0 - xe)*(x1 - xe)*(x2 - xe));
					m_cModH = (x2*x2*(h0*x1*x1*x1_mi_x2 - h1*x0*x0*x0_mi_x2) - he*x0_mi_x1*x0_mi_x2*x1_mi_x2*(x1*x2 + x0*(x1 + x2)) + (h1*(x0*x0*x0 - x2*x2*x2) + h0*(x2*x2*x2 - x1*x1*x1))*xe*xe + (h0*x1*x1 - h1*x0*x0 + (h1 - h0)*x2*x2)*xe*xe*xe + h2*x0_mi_x1*x0_mi_xe*x1_mi_xe*(x1*xe + x0*(x1 + xe)))*CommonMult6; // /((x0 - x1)*(x0 - x2)*(x1 - x2)*(x0 - xe)*(x1 - xe)*(x2 - xe));
					m_dModH = (he*x0*x0_mi_x1*x1*x0_mi_x2*x1_mi_x2*x2 + (-(h2*x0*x0_mi_x1*x1*x0_mi_xe*x1_mi_xe) + x2*(h1*x0*x0_mi_x2*x0_mi_xe - h0*x1*x1_mi_x2*x1_mi_xe)*x2_mi_xe)*xe)*CommonMult6; // /((x0 - x1)*(x0 - x2)*(x1 - x2)*(x0 - xe)*(x1 - xe)*(x2 - xe));
				}
				else
				{
					double CommonMult3 = inv_x0_mi_x1*inv_x0_mi_xe*inv_x1_mi_xe;
					m_bModH = (he*x0_mi_x1 + h0*x1_mi_xe - h1*x0_mi_xe)*CommonMult3;
					m_cModH = (he*(x1*x1 - x0*x0) + h1*(x0*x0 - xe*xe) + h0*(xe*xe - x1*x1))*CommonMult3;
					m_dModH = (he*x0*x0_mi_x1*x1 + xe*(h0*x1*x1_mi_xe - h1*x0*x0_mi_xe))*CommonMult3;
				}
			}
			else
			{
				m_cModH = (h0 - he)*inv_x0_mi_xe;
				m_dModH = (he*x0 - h0*xe)*inv_x0_mi_xe;
			}
		}
	}
}

//*************************************************************************
