/*
 * ProgramSettings.cpp
 *
 *  Created on: 18 бер. 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"

using namespace NonlinearFit;

Range readRange(const libconfig::Setting& stg)
{
	Range range;

	range.m_min = stg[0][0];
	range.m_max = stg[0][1];
	range.m_sampling = stg[1];

	return range;
}

FitParameter readFParameter(const libconfig::Setting& stg)
{
	bool isFit;
	static FitParameter fparameter;

	if(stg.isScalar())
	{
		fparameter.m_Name = stg.getPath();
		fparameter.m_Value = stg;
		fparameter.m_Lbvalue = 0.0;
		fparameter.m_Ubvalue = 0.0;
	}else if(stg.isGroup())
	{
		isFit = stg["tofit"];
		if(isFit)
		{
			fparameter.m_Name = stg.getPath();
			fparameter.m_Value = stg["value"];
			fparameter.m_Lbvalue = stg["lbvalue"];
			fparameter.m_Ubvalue = stg["ubvalue"];

			if(!stg.lookupValue("scvalue", fparameter.m_Scvalue))
			{
				fparameter.m_Scvalue = 0.0;
			}

			if(fparameter.m_Lbvalue >= fparameter.m_Ubvalue)
			{
				throw ProgramSettings::Exception("Fault boundaries " + toString(stg.getPath()));
			}
		}
		else
		{
			fparameter.m_Name = stg.getPath();
			fparameter.m_Value = stg["value"];
			fparameter.m_Lbvalue = 0.0;
			fparameter.m_Ubvalue = 0.0;
		}
	}
	else
	{
		throw ProgramSettings::Exception("Inapropriate type " + toString(stg.getPath()));
	}
	return fparameter;
}

std::ostream& operator<<(std::ostream& out, const Range& range)
{
	out << "[" << range.m_min << ", " << range.m_max << "]:" << range.m_sampling;
	return out;
}

std::ostream& operator<<(std::ostream& out, const FitParameter& fparam)
{
	/*variable parameter*/
	if((fparam.m_Lbvalue == 0.0) && (fparam.m_Ubvalue == 0.0))
	{
		out << fparam.m_Value;
	}
	else /*fixed parameter*/
	{
		out << fparam.m_Value << " [" << fparam.m_Lbvalue << " - " << fparam.m_Ubvalue << "]";
	}
	return out;
}

ProgramSettings::ProgramSettings()
{
}

ProgramSettings::~ProgramSettings()
{
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;

	m_cfgfile = cfgfile;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	/*reflection*/
	if(calculator["Q"].isArray() && calculator["Q"].getLength() == CalculatorSettings::HEXDIM)
	{
		for(int i = 0; i < CalculatorSettings::HEXDIM; ++i)
		{
			m_calculatorSettings.Q[i] = calculator["Q"][i];
		}
	}
	else
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}
	/*check the property of hexagonal Miller indices*/
	if((m_calculatorSettings.Q[0] + m_calculatorSettings.Q[1] + m_calculatorSettings.Q[2]) != 0)
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}

	/*X-ray wavelength*/
	m_calculatorSettings.lambda = calculator["lambda"];

	/*potential fit parameters*/
	m_calculatorSettings.scale = readFParameter(calculator["scale"]);
	m_calculatorSettings.background = readFParameter(calculator["background"]);
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	/*lattice parameters*/
	m_sampleSettings.a0 = sample["a0"];
	m_sampleSettings.c0 = sample["c0"];

	/*Poisson ratio*/
	m_sampleSettings.nu = sample["nu"];

	/*Sample thickness*/
	m_sampleSettings.thickness = sample["thickness"];

	/*dislocation settings*/
	const libconfig::Setting &dislocations = sample["dislocations"];
	m_sampleSettings.b_edge = dislocations["b_edge"];
	m_sampleSettings.b_screw = dislocations["b_screw"];

	/*potential fit parameters*/
	if(dislocations.exists("edge"))
	{
		m_sampleSettings.edgeDislocations.rho = readFParameter(dislocations["edge"]["rho"]);
		m_sampleSettings.edgeDislocations.Rc = readFParameter(dislocations["edge"]["rc"]);
	}
	else
	{
		m_sampleSettings.edgeDislocations.rho = FitParameter(dislocations.getPath() + ".edge.rho");
		m_sampleSettings.edgeDislocations.Rc = FitParameter(dislocations.getPath() + ".edge.rc");
	}

	if(dislocations.exists("screw"))
	{
		m_sampleSettings.screwDislocations.rho = readFParameter(dislocations["screw"]["rho"]);
		m_sampleSettings.screwDislocations.Rc = readFParameter(dislocations["screw"]["rc"]);
	}
	else
	{
		m_sampleSettings.screwDislocations.rho = FitParameter(dislocations.getPath() + ".screw.rho");
		m_sampleSettings.screwDislocations.Rc = FitParameter(dislocations.getPath() + ".screw.rc");
	}

	if(dislocations.exists("mixed"))
	{
		m_sampleSettings.mixedDislocations.rho = readFParameter(dislocations["mixed"]["rho"]);
		m_sampleSettings.mixedDislocations.Rc = readFParameter(dislocations["mixed"]["rc"]);
	}
	else
	{
		m_sampleSettings.mixedDislocations.rho = FitParameter(dislocations.getPath() + ".mixed.rho");
		m_sampleSettings.mixedDislocations.Rc = FitParameter(dislocations.getPath() + ".mixed.rc");
	}
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	defineMode(engine["mode"]);
	defineDiffractometry(engine["diffractometry"]);

	switch(m_engineSettings.workMode)
	{
	case EngineSettings::modeFITI:
		readIFitParameters(engine["parameters"]);
		break;
	case EngineSettings::modeCALCI:
	default:
		readICalcParameters(engine["parameters"]);
		break;
	}
	m_engineSettings.outfile = engine["outfile"].c_str();
}

void ProgramSettings::defineMode(const libconfig::Setting& stg)
{
	std::string mode;

	mode = stg.c_str();
	if (mode.compare("ICALC") == 0)
	{
		m_engineSettings.workMode = EngineSettings::modeCALCI;
	}
	else if(mode.compare("IFIT")==0)
	{
		m_engineSettings.workMode = EngineSettings::modeFITI;
	}
	else
	{
		m_engineSettings.workMode = EngineSettings::modeCALCI;
	}
}

void ProgramSettings::defineDiffractometry(const libconfig::Setting& stg)
{
	std::string diffractometry;

	diffractometry = stg.c_str();

	if (diffractometry.compare("DOUBLE") == 0)
	{
		m_engineSettings.workDiffractometry = EngineSettings::diffDOUBLE;
	}
	else if (diffractometry.compare("TRIPLE") == 0)
	{
		m_engineSettings.workDiffractometry = EngineSettings::diffTRIPLE;
	}
	else
	{
		m_engineSettings.workDiffractometry = EngineSettings::diffDOUBLE;
	}
}

void ProgramSettings::defineFitMethod(const libconfig::Setting& stg)
{
	std::string fit;

	fit = stg.c_str();
	if(fit.compare("LEVMAR") == 0)
	{
		m_engineSettings.fitMethod = EngineSettings::fitLEVMAR;
	}
	else if(fit.compare("PORT") == 0)
	{
		m_engineSettings.fitMethod = EngineSettings::fitPORT;
	}
	else if(fit.compare("NL2SOL") == 0)
	{
		/*TODO change when implemented*/
		m_engineSettings.fitMethod = EngineSettings::fitPORT;
	}
	else
	{
		m_engineSettings.fitMethod = EngineSettings::fitPORT;
	}
}

void ProgramSettings::defineFitType(const libconfig::Setting& stg)
{
	std::string fit;

	fit = stg.c_str();
	if(fit.compare("LOG") == 0)
	{
		m_engineSettings.fitType = NonlinearFitter::fitLOG;
	}
	else if(fit.compare("LIN") == 0)
	{
		m_engineSettings.fitType = NonlinearFitter::fitLIN;
	}
	else
	{
		m_engineSettings.fitType = NonlinearFitter::fitLIN;
	}
}

void ProgramSettings::readICalcParameters(const libconfig::Setting& stg)
{
	/*range of omega*/
	m_engineSettings.omegaRange = readRange(stg["omrange"]);
}

void ProgramSettings::readIFitParameters(const libconfig::Setting& stg)
{
	defineFitMethod(stg["method"]);
	defineFitType(stg["type"]);

	m_engineSettings.nbInterations = stg["nbIterations"];
	m_engineSettings.resumefile = stg["resumefile"].c_str();
	m_engineSettings.datafile = stg["datafile"].c_str();
}

void ProgramSettings::printEngineSettings() const
{
	std::cout << "---Engine settings---" << std::endl;

	switch(m_engineSettings.workMode)
	{
	case EngineSettings::modeCALCI:
		std::cout << "Mode: intensity calculation" << std::endl;
		printICalcParameters();
		break;
	case EngineSettings::modeFITI:
		std::cout << "Mode: intensity fit" << std::endl;
		printIFitParameters();
		break;
	default:
		//never happens
		break;
	}

	switch(m_engineSettings.workDiffractometry)
	{
	case EngineSettings::diffDOUBLE:
		std::cout << "Diffractometry: DOUBLE" << std::endl;
		break;
	case EngineSettings::diffTRIPLE:
		std::cout << "Diffractometry: TRIPLE" << std::endl;
		break;
	default:
		//never happens
		break;
	}

	std::cout << "Output basename:\t" << m_engineSettings.outfile << std::endl;
}

void ProgramSettings::printICalcParameters() const
{
	std::cout << "Omega range:\t" << m_engineSettings.omegaRange << std::endl;
}

void ProgramSettings::printIFitParameters() const
{
	switch (m_engineSettings.fitMethod)
	{
	case EngineSettings::fitLEVMAR:
		std::cout
				<< "Fit is performed by a version Levenberg-Marquadt algorithm (Levmar by M. Lourakis)"
				<< std::endl;
		break;
	case EngineSettings::fitNL2SOL:
		//break;
	case EngineSettings::fitPORT:
	default:
		std::cout
				<< "Fit is performed by a version of nl2sol algorithm (PORT by D.Gay)"
				<< std::endl;
		break;
	}
	switch (m_engineSettings.fitType)
	{
	case NonlinearFitter::fitLOG:
		std::cout
				<< "Logarithmic scale is used for fit"
				<< std::endl;
		break;
	case NonlinearFitter::fitLIN:
	default:
		std::cout
				<< "Linear scale is used for fit"
				<< std::endl;
		break;
	}
	std::cout << "Nb iterations:\t" << m_engineSettings.nbInterations
			<< std::endl;
	std::cout << "Data file:\t" << m_engineSettings.datafile << std::endl;
	std::cout << "Resume file:\t" << m_engineSettings.resumefile << std::endl;
}

void ProgramSettings::printSampleSettings() const
{
	std::cout << "---Sample settings---" << std::endl;
	std::cout << "Lattice parameters: (a0, c0)\t" << m_sampleSettings.a0 << ", "
			<< m_sampleSettings.c0 << std::endl;
	std::cout << "Sample thickness:\t" << m_sampleSettings.thickness << std::endl;
	std::cout << "Poisson ratio:\t" << m_sampleSettings.nu << std::endl;

	std::cout << "Edge Burgers vector:\t" << m_sampleSettings.b_edge << std::endl;
	std::cout << "Screw Burgers vector:\t" << m_sampleSettings.b_screw << std::endl;

	std::cout << "Screw dislocations:" << std::endl;
	std::cout << "\t-Density:\t" << m_sampleSettings.screwDislocations.rho << std::endl;
	std::cout << "\t-Correlation radius:\t" << m_sampleSettings.screwDislocations.Rc << std::endl;

	std::cout << "Edge dislocations:\t" << std::endl;
	std::cout << "\t-Density:\t" << m_sampleSettings.edgeDislocations.rho << std::endl;
	std::cout << "\t-Correlation radius:\t" << m_sampleSettings.edgeDislocations.Rc << std::endl;

	std::cout << "Mixed dislocations:\t" << std::endl;
	std::cout << "\t-Density:\t" << m_sampleSettings.mixedDislocations.rho << std::endl;
	std::cout << "\t-Correlation radius:\t" << m_sampleSettings.mixedDislocations.Rc << std::endl;
}

void ProgramSettings::printCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;
	std::cout << "Reflection:\t[" << m_calculatorSettings.Q[0] << ", "
			<< m_calculatorSettings.Q[1] << ", " << m_calculatorSettings.Q[2]
			<< ", " << m_calculatorSettings.Q[3] << "]" << std::endl;
	std::cout << "X-ray wavelength:\t" << m_calculatorSettings.lambda << std::endl;
	std::cout << "Intensity scale coefficient:\t" << m_calculatorSettings.scale << std::endl;
	std::cout << "Intensity background:\t" << m_calculatorSettings.background << std::endl;
}

void ProgramSettings::print() const
{
	printEngineSettings();
	printSampleSettings();
	printCalculatorSettings();
}
