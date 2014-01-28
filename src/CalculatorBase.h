/*
 * CalculatorBase.h
 *
 *  Created on: 15 june 2012
 *      Author: kopp
 */

#ifndef CALCULATORBASE_H_
#define CALCULATORBASE_H_

#include "NonlinearFit.h"
#include <StringTools.h>
#include <Log.h>

#include <vector>

class CalculatorBase: public NonlinearFit::FitCalculator
{
public:
	CalculatorBase(double Qx, double Qz, double be, double bs, double lambda, double nu, double l);
	virtual ~CalculatorBase();

	virtual void init(const double * stgs) {};
	virtual double getI(const double& q) = 0;
	virtual double getG(const double& r);

	virtual void reinit(const NonlinearFit::CalculatorParameterMap& params);
	double eval(const NonlinearFit::CalculatorArgument * arg);
protected:

	double L;//sample thickness
	double m_Qnorm, m_Qx, m_Qz;

	double Isc, Ibg;
	double rho_edge, rho_screw, rho_mixed,
			rc_edge, rc_screw, rc_mixed;

	double C_edge, C_screw, C_mixed;
	double ksi_mixed, ksi_screw, ksi_edge;
	double cosPhi, tanPhi, cosThetaB, sinThetaB;
};

class CalculatorBaseArgument : public NonlinearFit::CalculatorArgument
{
public:
	CalculatorBaseArgument(double val = 0) {m_Value = val;}
	virtual ~CalculatorBaseArgument() {}
	double m_Value;
};

#endif /* CALCULATORBASE_H_ */
