/*
 * LevmarFitter.h
 *
 *  Created on: 29 квіт. 2013
 *      Author: kopp
 */

#ifndef LEVMARFITTER_H_
#define LEVMARFITTER_H_

#include <levmar.h>
#include "NonlinearFit.h"

namespace NonlinearFit
{

class LevmarFitter: public NonlinearFit::NonlinearFitter
{
public:
	LevmarFitter();
	virtual ~LevmarFitter();
	virtual void fit(FitType type = fitLIN, int nbit=100);
	virtual std::string getReasonToStop() const;
protected:
	friend void LinNormFuncLM(double *x, double *v, int m, int n, void * adata);
	friend void LogNormFuncLM(double *x, double *v, int m, int n, void * adata);

	void resetCalculator(const double * x);
	ResidualValueType getCalcResidual(size_t i);

	void extractFitInfo(const double * info, const double * covar);
};

} /* namespace NonlinearFit */
#endif /* LEVMARFITTER_H_ */
