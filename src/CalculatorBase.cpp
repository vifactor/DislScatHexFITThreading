/*
 * CalculatorBase.cpp
 *
 *  Created on: 15 june 2012
 *      Author: kopp
 */

#include "CalculatorBase.h"

using namespace NonlinearFit;

CalculatorBase::CalculatorBase(double Qx, double Qz, double be, double bs, double lambda, double nu, double l) :
		FitCalculator()
{
	L = l;
	m_Qx = Qx;
	m_Qz = Qz;
	m_Qnorm = sqrt(m_Qx * m_Qx + m_Qz * m_Qz);

	/* Bragg's law */
	sinThetaB = m_Qnorm * lambda / (4 * M_PI);
	cosThetaB = sqrt(1 - gsl_pow_2(sinThetaB));

	/* Rc = L / ksi, where ksi = (Q.b)/2pi*/
	if(m_Qz !=0 )
		ksi_screw = fabs(m_Qz * bs / (2 * M_PI));
	else
		ksi_screw = 1.0;

	if(m_Qx != 0)
		ksi_edge = fabs(m_Qx * be / (2 * M_PI));
	else
		ksi_edge = 1.0;
	ksi_mixed = (m_Qz * bs + m_Qx * be) / (2 * M_PI);

	C_screw = 0.0;
	C_edge = 0.0;
	C_mixed = 0.0;
	rho_edge = 0.0;
	rho_screw = 0.0;
	rho_mixed = 0.0;
	rc_edge = 0.0;
	rc_screw = 0.0;
	rc_mixed = 0.0;
	cosPhi = 0.0;
	tanPhi = 0.0;
	Isc = 1.0;
	Ibg = 0.0;

	std::cout << "Qnorm:\t" << m_Qnorm << std::endl;
	std::cout << "sin(thetaB):\t" << sinThetaB << std::endl;
	std::cout << "ThetaB [deg]:\t" << asin(sinThetaB) * 180 / M_PI << std::endl;
}

CalculatorBase::~CalculatorBase()
{
}

double CalculatorBase::getG(const double& r)
{
	static double T_edge, T_screw, T_mixed;

	T_edge = 0;
	if((rho_edge > 0) && (fabs(C_edge) > 0))
	{
		T_edge = C_edge * rho_edge * r * r *
				log((rc_edge / ksi_edge + r) / r);
	}

	T_screw = 0;
	if((rho_screw > 0) && (fabs(C_screw) > 0))
	{
		T_screw = C_screw * rho_screw * r * r *
				log((rc_screw / ksi_screw + r) / r);
	}

	T_mixed = 0;
	if(rho_mixed > 0)
	{
		T_mixed = C_mixed * rho_mixed * r * r *
				log((rc_mixed / ksi_mixed + r) / r);
	}

	//std::cout << "r:\t" << r << std::endl;
	//std::cout << "C_mixed * rho_mixed:\t" << C_mixed * rho_mixed << std::endl;
	//std::cout << "C_edge * rho_edge:\t" << C_edge * rho_edge << std::endl;
	//std::cout << "C_screw * rho_screw:\t" << C_screw * rho_screw << std::endl;
	//std::cout << "ksi_screw:\t" << ksi_screw << std::endl;
	//std::cout << "ksi_edge:\t" << ksi_edge << std::endl;
	//std::cout << "ksi_mixed:\t" << ksi_mixed << std::endl;
	//std::cout << "Tscrew:\t" << T_screw << std::endl;
	//std::cout << "Tedge:\t" << T_edge << std::endl;
	//std::cout << "Tmixed:\t" << T_mixed << std::endl;

	return exp(-(T_edge + T_screw + T_mixed));
}


void CalculatorBase::reinit(const NonlinearFit::CalculatorParameterMap& params)
{
	Isc = params.find("Calculator.scale")->second;
	Ibg = params.find("Calculator.background")->second;

	/*initially dislocation density is given in [cm-2]*/
	rho_edge = params.find("Sample.dislocations.edge.rho")->second  * 1e-14;
	rho_screw = params.find("Sample.dislocations.screw.rho")->second * 1e-14;
	rho_mixed = params.find("Sample.dislocations.mixed.rho")->second * 1e-14;

	rc_edge = params.find("Sample.dislocations.edge.rc")->second;
	rc_screw = params.find("Sample.dislocations.screw.rc")->second;
	rc_mixed = params.find("Sample.dislocations.mixed.rc")->second;

	//std::cout << "rho_edge:\t" << rho_edge << std::endl;
	//std::cout << "rc_edge:\t" << rc_edge << std::endl;
	//std::cout << "Isc:\t" << Isc << std::endl;
}

double CalculatorBase::eval(const NonlinearFit::CalculatorArgument * arg)
{
	return getI(static_cast<const CalculatorBaseArgument* >(arg)->m_Value);
}
