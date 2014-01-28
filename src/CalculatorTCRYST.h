/*
 * CalculatorTCRYST.h
 *
 *  Created on: 10 ρεπο. 2012
 *      Author: kopp
 */

#ifndef CALCULATORTCRYST_H_
#define CALCULATORTCRYST_H_

#include "CalculatorBase.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

class CalculatorTCRYST: public CalculatorBase
{
public:
	enum {stgsNB = 1};
	CalculatorTCRYST(double Qx, double Qz, double be,	double bs, double lambda, double nu, double l);
	double getI(const double& q);
	virtual ~CalculatorTCRYST();
protected:
	gsl_integration_workspace * ws;
	gsl_function F;
	size_t limit;
	friend double tcryst_integrand_asym(double var, void *params);
	friend double tcryst_integrand_sym(double var, void *params);

	double q;
};

#endif /* CALCULATORTCRYST_H_ */
