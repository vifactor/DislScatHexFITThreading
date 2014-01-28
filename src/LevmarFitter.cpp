/*
 * LevmarFitter.cpp
 *
 *  Created on: 29 квіт. 2013
 *      Author: kopp
 */

#include "LevmarFitter.h"

namespace NonlinearFit
{

LevmarFitter::LevmarFitter() : NonlinearFitter()
{
	m_Name = "LEVMAR";
}

LevmarFitter::~LevmarFitter()
{
}

void LinNormFuncLM(double * x, double * v, int m, int n, void * adata)
{
	/* adata is a pointer to the fitter */
	static LevmarFitter * fitter;

	fitter = reinterpret_cast<LevmarFitter *>(adata);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */
	for (int i = 0; i < n; ++i)
	{
		v[i] = fitter->getCalcResidual(i);
	}
}

void LogNormFuncLM(double * x, double * v, int m, int n, void * adata)
{
	register int i;
	register double resid;

	/* adata is a pointer to the fitter */
	LevmarFitter * fitter = reinterpret_cast<LevmarFitter *>(adata);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */
	for (i = 0; i < n; ++i)
	{
		resid = fitter->getCalcResidual(i);
		if(resid > 0)
		{
			v[i] = log(resid);
		}
		else
		{
			v[i] = 0.0;
		}
	}
}

void LevmarFitter::fit(FitType type, int nbit)
{
	int nbParams, nbResids;
	/* levmar specific parameters */
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	double * x_ptr;
	double * xlb_ptr;
	double * xub_ptr;
	double * xdscl_ptr;
	double * covar_ptr;
	double * residuals;
	void (*NormFunc)(double *x, double *v, int m, int n, void * adata);

	nbParams = m_FitParameters.size();
	nbResids = m_DataPoints.size();

	/*allocate arrays*/
	x_ptr = new double[nbParams];
	xlb_ptr = new double[nbParams];
	xub_ptr = new double[nbParams];
	xdscl_ptr = new double[nbParams];
	covar_ptr = new double[nbParams * nbParams];
	residuals = new double[nbResids];

	std::cout << "Allocation" << std::endl;

	/*select norm function pointer*/
	if(type == fitLIN)
	{
		NormFunc = LinNormFuncLM;
	}else if(type == fitLOG)
	{
		NormFunc = LogNormFuncLM;
		for(size_t i = 0; i < m_DataPoints.size(); ++i)
			m_DataPoints[i].m_Residual = log(m_DataPoints[i].m_Residual);
	}
	else
	{
		NormFunc = LinNormFuncLM;
	}

	//initialize levmar calc options
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-1;
	opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

	//initialize levmar parameters
	for(size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		x_ptr[i]	= m_FitParameters[i].m_Value;
		xlb_ptr[i]	= m_FitParameters[i].m_Lbvalue;
		xub_ptr[i]	= m_FitParameters[i].m_Ubvalue;
		xdscl_ptr[i]= m_FitParameters[i].m_Scvalue;
	}

	/*setup residuals*/
	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		residuals[i] = m_DataPoints[i].m_Residual;
	}

	dlevmar_bc_dif(NormFunc, x_ptr, residuals, nbParams,
			nbResids, xlb_ptr, xub_ptr, xdscl_ptr, nbit,
			opts, info, NULL, covar_ptr, reinterpret_cast<void *>(this));

	resetFitParameters(x_ptr);
	extractFitInfo(info, covar_ptr);

	/*deallocate arrays*/
	delete[] x_ptr;
	delete[] xlb_ptr;
	delete[] xub_ptr;
	delete[] xdscl_ptr;
	delete[] covar_ptr;
	delete[] residuals;
}

void LevmarFitter::resetCalculator(const double * x)
{
	NonlinearFitter::resetCalculatorParameters(x);
	m_Calculator->reinit(m_CalculatorParameters);
}

ResidualValueType LevmarFitter::getCalcResidual(size_t i)
{
	return m_Calculator->eval(m_DataPoints[i].m_Argument);
}

void LevmarFitter::extractFitInfo(const double * info, const double * covar)
{
	size_t index;

	m_Finit = sqrt(info[0]);
	m_Ffin = sqrt(info[1]);
	m_nbFuncEval = info[7];
	m_nbGradEval = info[8];
	m_nbIterPerf = info[5];
	m_reasonToStopId = info[6];

	/*extract covariance matrix*/
	index = 0;
	for(size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		for(size_t j = 0; j <= i; ++j)
		{
			index = j + i * m_FitParameters.size();
			m_matrixCovar[i][j] = covar[index];
		}
	}
}

std::string LevmarFitter::getReasonToStop() const
{
	switch(m_reasonToStopId)
	{
	case 1:
		return "Small gradient J^T f";
		break;
	case 2:
		return "Small Dp";
		break;
	case 3:
		return "Max nb iterations";
		break;
	case 4:
		return "Singular matrix";
		break;
	case 5:
		return "No further error reduction is possible";
		break;
	case 6:
		return "Small ||f||^2";
		break;
	default:
		return "Invalid parameters";
		break;
	}
	return "Invalid parameters";
}

} /* namespace NonlinearFit */
