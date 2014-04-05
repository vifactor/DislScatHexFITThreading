/*
 * CalculatorTCRYST.cpp
 *
 *  Created on: 10 august 2012
 *      Author: kopp
 */

#include "CalculatorTCRYST.h"

double tcryst_integrand_sym(double var, void *params)
{
	static double result = 0;
	CalculatorTCRYST * me = static_cast<CalculatorTCRYST *>(params);

	result = me->getG(var) * cos(me->q * var);

	//std::cout << "G(" << var << ")=" << result;
	return result;
}

double tcryst_integrand_asym(double var, void *params)
{
	static double result = 0;
	CalculatorTCRYST * me = static_cast<CalculatorTCRYST *>(params);

	result = me->getG(var) * var * gsl_sf_bessel_J0(me->q * var);

	return result;
}

CalculatorTCRYST::CalculatorTCRYST(double Qx, double Qz, double be, double bs,
		double lambda, double nu, double l) :
		CalculatorBase(Qx, Qz, be, bs, lambda, nu, l)
{
	double cosPsi, sinPsi, sinPhi;
	double gamma_edge, gamma_screw;

	q = 0.0;

	/* Psi - angle between (x, y) plane and Q
	 *
	 * (Q.n) = |Q| |n| cos(Pi/2 - Psi), n = (0,0,1) => sinPsi = Qz / Q
	 */
	sinPsi = m_Qz / m_Qnorm;
	cosPsi = sqrt(1 - gsl_pow_2(sinPsi));

	/* Phi  - angle between (x, y) plane and Kout
	 *
	 * can be estimated from Bragg's law and relation (Kout - Kin) = Q
	 *
	 */
	sinPhi = sinPsi * sinThetaB;
	cosPhi = sqrt(1 - gsl_pow_2(sinPhi));

	std::cout << "sin(Psi):\t" << sinPsi << std::endl;
	std::cout << "cos(Phi):\t" << cosPhi << std::endl;
	std::cout << "sin(Phi):\t" << sinPhi << std::endl;

	/*average value of chi_edge*/
	gamma_edge = (3 - 6 * nu + 4 * nu * nu)
			* gsl_pow_2(cosPsi) / (8 * gsl_pow_2(nu - 1));

	/*average value of chi_screw*/
	gamma_screw = sinPsi * sinPsi / 2;

	C_edge = gamma_edge * gsl_pow_2(m_Qnorm * be) / (4 * M_PI);
	C_screw = gamma_screw * gsl_pow_2(m_Qnorm * bs) / (4 * M_PI);
	C_mixed = C_edge + C_screw;


	std::cout << "C_edge:\t" << C_edge << std::endl;
	std::cout << "C_screw:\t" <<C_screw << std::endl;
	std::cout << "C_mixed:\t" <<C_mixed << std::endl;

	/*--------integration stuff------*/
	limit = 20000;
	F.params = this;
	if(sinPsi == 1)
	{
		std::cout << "Symmetric reflection" << std::endl;
		F.function = &tcryst_integrand_sym;
	}
	else
	{
		std::cout << "Asymmetric reflection" << std::endl;
		F.function = &tcryst_integrand_asym;
	}

	ws = gsl_integration_workspace_alloc(limit);

	gsl_set_error_handler_off ();
}

CalculatorTCRYST::~CalculatorTCRYST()
{
	if(ws)
		gsl_integration_workspace_free(ws);
}

double CalculatorTCRYST::getI(const double& om)
{
	static double result;
	static double error;
	static const double epsabs = 0.0, epsrel = 1e-6;

	/*
	 *
	 * q = Q0 * omega
	 *
	 */
	q = om * m_Qnorm;
	result = 0.0;
	error = 0.0;
	gsl_integration_qagiu (&F, 0.0, epsabs, epsrel, limit, ws, &result, &error);

	//LOG(logDEBUG1) << "I(" << q << ")=" << result;

	return Isc * result + Ibg;
}
