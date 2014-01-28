/*
 * ProgramSettings.h
 *
 *  Created on: 18 бер. 2013
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include <exception>
#include <libconfig.h++>
#include <StringTools.h>
#include "UsefulStructures.h"
#include "NonlinearFit.h"

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "ProgramSettings::" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	struct SampleSettings
	{
		double nu;
		double thickness;
		/*hexagonal lattice parameters*/
		double a0, c0;
		/*burgers components*/
		double b_edge, b_screw;

		struct
		{
			/*density*/
			NonlinearFit::FitParameter rho;
			/*correlation radius*/
			NonlinearFit::FitParameter Rc;
		}
		edgeDislocations, screwDislocations, mixedDislocations;
	};
	struct CalculatorSettings
	{
		enum {HEXDIM = 4};
		/*hexagonal index*/
		int Q[HEXDIM];
		/*wavelength*/
		double lambda;

		NonlinearFit::FitParameter scale;
		NonlinearFit::FitParameter background;
	};
	struct EngineSettings
	{
		std::string outfile;
		enum WorkDiffractometry {diffDOUBLE, diffTRIPLE} workDiffractometry;
		enum WorkMode {modeCALCI, modeFITI}	workMode;

		/*only for intensity calculation*/
		Range omegaRange;

		/*only for intensity fit*/
		enum FitMethod {fitLEVMAR, fitPORT, fitNL2SOL} fitMethod;
		NonlinearFit::NonlinearFitter::FitType fitType;
		size_t nbInterations;
		std::string datafile;
		std::string resumefile;
	};
	ProgramSettings();
	virtual ~ProgramSettings();
	const SampleSettings& getSampleSettings() const
	{
		return m_sampleSettings;
	}
	const CalculatorSettings& getCalculatorSettings() const
	{
		return m_calculatorSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return m_engineSettings;
	}
	const std::string& getConfigfile() const
	{
		return m_cfgfile;
	}
	void read(const std::string& cfg);
	void print() const;
	protected:
	void readCalculatorSettings(const libconfig::Setting& root);
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void defineMode(const libconfig::Setting& stg);
	void defineDiffractometry(const libconfig::Setting& stg);
	void defineFitMethod(const libconfig::Setting& stg);
	void defineFitType(const libconfig::Setting& stg);

	void readIFitParameters(const libconfig::Setting& stg);
	void readICalcParameters(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printCalculatorSettings() const;
	void printEngineSettings() const;

	void printIFitParameters() const;
	void printICalcParameters() const;

	SampleSettings m_sampleSettings;
	CalculatorSettings m_calculatorSettings;
	EngineSettings m_engineSettings;
	std::string m_cfgfile;
};

#endif /* PROGRAMSETTINGS_H_ */
