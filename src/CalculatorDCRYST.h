/*
 * CalculatorSQUAD.h
 *
 *  Created on: 28 june 2012
 *      Author: kopp
 */

#ifndef CALCULATORDCRYST_H_
#define CALCULATORDCRYST_H_

#include "CalculatorBase.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>

class CalculatorDCRYST : public CalculatorBase
{
public:
	enum {stgsNBINTRP, stgsMINOM, stgsMAXOM, stgsNB};
	CalculatorDCRYST(double Qx, double Qz, double be, double bs,
			double lambda, double nu, double l);
	double getI(const double& q);
	virtual double getG(const double& r);
	virtual ~CalculatorDCRYST();
private:
	gsl_integration_workspace * ws;
	gsl_integration_workspace * wsc;
	gsl_integration_qawo_table * qawo_table;
	gsl_function F, F_qq;
	size_t limit;
	friend double dcryst_integrand(double var, void *params);
	friend double dcryst_integrand_qq(double var, void *params);

	double qq;
};

#endif /* CALCULATORDCRYST_H_ */
