/*
 * CalculatorDCRYST.cpp
 *
 *  Created on: 28 june 2012
 *      Author: kopp
 */

#include "CalculatorDCRYST.h"

double dcryst_integrand(double var, void *params)
{
	static double result = 0;

	result = static_cast<CalculatorDCRYST *>(params)->getG(var);
	return result;
}

double dcryst_integrand_qq(double var, void *params)
{
	static double result = 0;

	CalculatorDCRYST * me = static_cast<CalculatorDCRYST *>(params);

	result = me->getG(var) * cos(var * me->qq);
	return result;
}

CalculatorDCRYST::CalculatorDCRYST(double Qx, double Qz, double be, double bs,
		double lambda, double nu, double l) :
		CalculatorBase(Qx, Qz, be, bs, lambda, nu, l)
{
	double cosPsi, sinPsi, cosAlpha, sinPhi;
	double gamma_edge, gamma_screw;

	qq = 0.0;
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
	tanPhi = sinPhi / cosPhi;

	/*
	 * alpha - angle between projections of Q and Kout on (x, y) plane
	 *
	 */
	cosAlpha = sinThetaB * cosPsi / cosPhi;

	LOG(logINFO) << "cos(Alpha):\t" << cosAlpha << std::endl;
	LOG(logINFO) << "sin(Psi):\t" << sinPsi << std::endl;
	LOG(logINFO) << "cos(Phi):\t" << cosPhi << std::endl;

	/*precise value of chi_edge in case of double crystal diffractometry*/
	gamma_edge = (9 - 16 * nu + 8 * nu * nu - 2 * (3 - 4 * nu) * gsl_pow_2(cosAlpha)) * gsl_pow_2(cosPsi)
				/ (16 * gsl_pow_2(nu - 1));

	LOG(logDEBUG) << "gamma_edge (correct):\t" << gamma_edge << std::endl;

	gamma_screw = gsl_pow_2(sinPsi)/ 2;

	C_edge = gamma_edge * gsl_pow_2(m_Qnorm * be) / (4 * M_PI);
	C_screw = gamma_screw * gsl_pow_2(m_Qnorm * bs) / (4 * M_PI);
	C_mixed = (C_edge + C_screw);

	LOG(logDEBUG) << "gamma_edge:\t" <<gamma_edge << std::endl;
	LOG(logDEBUG) << "gamma_screw:\t" <<gamma_screw << std::endl;

	LOG(logDEBUG) << "gamma_screw/gamma_edge:\t" <<gamma_edge / gamma_screw << std::endl;

	LOG(logDEBUG) << "C_edge:\t" <<C_edge << std::endl;
	LOG(logDEBUG) << "C_screw:\t" <<C_screw << std::endl;
	LOG(logDEBUG) << "C_mixed:\t" <<C_mixed << std::endl;

	/*------integration stuff--------*/
	limit = 20000;
	F.function 		= &dcryst_integrand;
	F.params 		= this;

	F_qq.function 		= &dcryst_integrand_qq;
	F_qq.params 		= this;

	qawo_table 		= gsl_integration_qawo_table_alloc(0.0, 0.0, GSL_INTEG_COSINE, limit);
	ws 				= gsl_integration_workspace_alloc(limit);
	wsc				= gsl_integration_workspace_alloc(limit);

	gsl_set_error_handler_off ();
}

CalculatorDCRYST::~CalculatorDCRYST()
{
	if(qawo_table)
		gsl_integration_qawo_table_free(qawo_table);
	if(wsc)
		gsl_integration_workspace_free(wsc);
	if(ws)
		gsl_integration_workspace_free(ws);
}

double CalculatorDCRYST::getI(const double& om)
{
	static double result;
	static double error;
	static const double epsabs = 0.0, epsrel = 1e-6;

	qq = om * m_Qnorm * cosThetaB / cosPhi;

	gsl_integration_qagiu (&F_qq, 0.0, epsabs, epsrel, limit, ws, &result, &error);

	//gsl_integration_qawo_table_set(qawo_table, qq, 0, GSL_INTEG_COSINE);
	//gsl_integration_qawf (&F, 0, epsabs, limit, ws, wsc, qawo_table, &iIntrp[i], &error);

	return Isc * result + Ibg;
}

double CalculatorDCRYST::getG(const double& r)
{
	static double L_correct;

	L_correct = 1.0;
	/*if (cosPhi != 0)
	{
		L_correct = exp(-2 * tanPhi / L * r);
	}*/

	return CalculatorBase::getG(r) * L_correct;
}
