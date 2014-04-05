/*
 *
 *  Created on: 30 may 2012
 *      Author: kopp
 */

#include "Engine.h"

using namespace NonlinearFit;

Engine::Engine()
{
	m_calculator = NULL;
	m_fitter = NULL;
	m_programSettings = NULL;
}

Engine::~Engine()
{
	if(m_calculator)
		delete m_calculator;
	if(m_fitter)
		delete m_fitter;
	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		if(m_DataPoints[i].m_Argument)
			delete m_DataPoints[i].m_Argument;
	}
}

void Engine::exec(const ProgramSettings * stg)
{
	m_programSettings = stg;

	m_programSettings->print();
	init();
	setupComponents();
	doWork();
	saveResult();
}

void Engine::init()
{
	if(m_calculator)
		delete m_calculator;
	if(m_fitter)
		delete m_fitter;
	/*vector of fit parameters with their initial values*/
	m_fParametersInit.clear();
	/*vector of fit parameters with their final values*/
	m_fParametersFinal.clear();
	/*map of calculator parameters (potential fit parameters)*/
	m_cParameters.clear();

	calcValX.clear(); calcValZini.clear(); calcValZfin.clear();
	origValX.clear(); origValZ.clear();

	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		delete m_DataPoints[i].m_Argument;
	}
	m_DataPoints.clear();
}

void Engine::doWork()
{
	//m_calculator->reinit(m_cParameters);
	//m_calculator->getG(1);

	switch (m_programSettings->getEngineSettings().workMode)
	{
	case ProgramSettings::EngineSettings::modeCALCI:
		calcI(calcValZini);
		break;
	case ProgramSettings::EngineSettings::modeFITI:
		calcI(calcValZini);
		fitI();
		printFitterInfo(m_fitter);
		calcI(calcValZfin);
		saveSettings();
		break;
	default:
		//this should never happen
		break;
	}
}

void Engine::calcI(std::vector<double>& vals)
{
	double q, I;

	m_calculator->reinit(m_cParameters);

	for(size_t ipt = 0; ipt < vals.size(); ipt++)
	{
		q = calcValX[ipt];
		I = m_calculator->getI(q);

		//LOG(logDEBUG) << q << "\t" << I << std::endl;

		vals[ipt] = I;
	}
}

std::string Engine::getFilename() const
{
	std::string filename;
	std::string extension;

	filename = m_programSettings->getEngineSettings().outfile;
	extension = stripExtension(filename);
	switch (m_programSettings->getEngineSettings().workMode)
	{
	case ProgramSettings::EngineSettings::modeCALCI:
		filename += "_i";
		break;
	case ProgramSettings::EngineSettings::modeFITI:
		filename += "_ft";
		break;
	default:
		break;
	}
	filename += "." + extension;
	return filename;
}

void Engine::saveResult() const
{
	std::ofstream fout(getFilename().c_str());
	if(!fout)
	{
		throw Engine::Exception("Unable to open the output file:\t" + getFilename());
	}

	fout<<"#domega[deg]\tdomega[rad]\tIexp\tIini\tIfin"<<std::endl;
	for(size_t ipt = 0; ipt < origValZ.size(); ipt++)
	{
		fout<<origValX[ipt]<<"\t"
			<<calcValX[ipt]<<"\t"
			<<origValZ[ipt]<<"\t"
			<<calcValZini[ipt]<<"\t"
			<<calcValZfin[ipt]<<"\t"
				<<std::endl;
	}
	fout.close();
	switch (m_programSettings->getEngineSettings().workMode)
	{
	case ProgramSettings::EngineSettings::modeCALCI:
		break;
	case ProgramSettings::EngineSettings::modeFITI:
		saveResume();
		break;
	default:
		break;
	}

}

void Engine::fitI()
{
	//carry out fit
	//with max "calculationSettings.modefitSettings.nbIterations" iterations
	m_fitter->fit(m_programSettings->getEngineSettings().fitType, m_programSettings->getEngineSettings().nbInterations);

	//get new values (best fit values) of fit parameters
	m_fParametersFinal = m_fitter->getFitParameters();

	//update calc parameters with new fitted values
	mergeParameters(m_cParameters, m_fParametersFinal);
}

void Engine::saveResume() const
{
	double rho_screw, rho_edge, rho_mixed,
		rc_screw, rc_edge, rc_mixed,
		M_edge, M_screw, M_mixed;
	std::string filename;
	std::ofstream fout, fin;

	rc_edge = m_cParameters.find("Sample.dislocations.edge.rc")->second;
	rc_screw = m_cParameters.find("Sample.dislocations.screw.rc")->second;
	rc_mixed = m_cParameters.find("Sample.dislocations.mixed.rc")->second;
	rho_edge = m_cParameters.find("Sample.dislocations.edge.rho")->second;
	rho_screw = m_cParameters.find("Sample.dislocations.screw.rho")->second;
	rho_mixed = m_cParameters.find("Sample.dislocations.mixed.rho")->second;

	/*transform rc to M = rc/rd = rc*rho^(1/2) */
	M_screw = rc_screw * sqrt(rho_screw * 1e-14);
	M_edge = rc_edge * sqrt(rho_edge * 1e-14);
	M_mixed = rc_mixed * sqrt(rho_mixed * 1e-14);

	filename = m_programSettings->getEngineSettings().resumefile;

	fin.open(filename.c_str(), std::ios::in);
	if(fin.good())
	{
		fin.close();
		fout.open(filename.c_str(), std::ios::app);
	}
	else
	{
		fin.close();
		fout.open(filename.c_str(), std::ios::out);
		fout << "#cfg\tdat\tH\tK\tI\tL\t"<<
				"rho_edge\trc_edge\tM_edge\t" <<
				"rho_screw\trc_screw\tM_screw\t" <<
				"rho_mixed\trc_mixed\tM_mixed\n";
	}
	fout << m_programSettings->getConfigfile() << "\t";
	fout << m_programSettings->getEngineSettings().datafile << "\t";
	fout << m_programSettings->getCalculatorSettings().Q[0] << "\t" <<
			m_programSettings->getCalculatorSettings().Q[1] << "\t" <<
			m_programSettings->getCalculatorSettings().Q[2] << "\t" <<
			m_programSettings->getCalculatorSettings().Q[3] << "\t";
	fout << rho_edge << "\t";
	fout << rc_edge << "\t";
	fout << M_edge << "\t";
	fout << rho_screw << "\t";
	fout << rc_screw << "\t";
	fout << M_screw << "\t";
	fout << rho_mixed << "\t";
	fout << rc_mixed << "\t";
	fout << M_mixed << "\t";
	fout << std::endl;
	fout.close();
}

void Engine::saveSettings() const
{
	libconfig::Config cfg;
	std::string oldcfgfile, newcfgfile;

	oldcfgfile = m_programSettings->getConfigfile();
	newcfgfile = m_programSettings->getConfigfile();
	stripExtension(newcfgfile);
	newcfgfile += "_mod.cfg";

	try
	{
		// Read the configuration file. If there is an error, report it
		cfg.readFile(oldcfgfile.c_str());

		//reset configuration putting the new fitted values from fitParameterListFinal
		//for(auto fparameter : fitParameterListFinal)
		for(size_t i = 0; i < m_fParametersFinal.size(); ++i)
		{
			cfg.lookup(m_fParametersFinal[i].m_Name + ".value") = m_fParametersFinal[i].m_Value;
		}
		//save new configuration
		cfg.writeFile (newcfgfile.c_str());
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Engine::Exception("I/O error while reading file:\t" + oldcfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Engine::Exception("Parse error at " +
									toString(pex.getFile()) + ":" +
									toString(pex.getLine()) + " - " +
									toString(pex.getError()));
	}catch(const libconfig::SettingNotFoundException &nfex)
	{
		throw Engine::Exception(toString(nfex.getPath()));
	}catch(libconfig::SettingTypeException& tex)
	{
		throw Engine::Exception(toString(tex.getPath()) + "(" + toString(tex.what()) + ")");
	}
}

void Engine::setupComponents()
{
	setupCalculator();
	setupCParameters();
	switch(m_programSettings->getEngineSettings().workMode)
	{
		case ProgramSettings::EngineSettings::modeCALCI:
			fillData();
			break;
		case ProgramSettings::EngineSettings::modeFITI:
			readData();
			setupFitter();
			break;
		default:
			//never happens
			break;
	}
}

void Engine::fillData()
{

	m_programSettings->getEngineSettings().omegaRange.toVector(origValX);
	m_programSettings->getEngineSettings().omegaRange.toVector(calcValX, M_PI / 180.0);

	origValZ.resize(origValX.size(), 0.0);
	calcValZini.resize(origValX.size(), 0.0);
	calcValZfin.resize(origValX.size(), 0.0);
}

void Engine::readData()
{
	DataReader dr;
	dr.readFile(m_programSettings->getEngineSettings().datafile);

	if(!dr.good())
	{
		throw Engine::Exception("File " + m_programSettings->getEngineSettings().datafile +
				" has not been read or has no a header");
	}

	if(dr.columnExist("[intensity]"))
	{
		dr.getColumn(origValZ, "[intensity]");
	}
	else
	{
		throw Engine::Exception("Column \"[intensity]\" has not been found in " +
				m_programSettings->getEngineSettings().datafile);
	}

	if (dr.columnExist("[domega][rad]"))
	{
		//get points without any transformation
		dr.getColumn(origValX, "[domega][rad]");
		/*allocate space*/
		calcValX.resize(origValX.size());
		/*copy arrays*/
		std::copy(origValX.begin(), origValX.end(), calcValX.begin());
	}
	else if (dr.columnExist("[omega][rad]"))
	{
		//get points without any transformation
		dr.getColumn(origValX, "[omega][rad]");
		/*allocate space*/
		calcValX.resize(origValX.size());
		/*copy to radians*/
		std::copy(origValX.begin(), origValX.end(), calcValX.begin());
		/*subtracts the x_max value from array*/
		toCenter(calcValX, origValZ);
	}
	else if (dr.columnExist("[domega]"))
	{
		//get points without any transformation
		dr.getColumn(origValX, "[domega]");
		/*allocate space*/
		calcValX.resize(origValX.size());
		/*transform to radians*/
		transform(origValX.begin(), origValX.end(), calcValX.begin(),
		          bind2nd(std::multiplies<double>(), M_PI / 180.0));
	}
	else if (dr.columnExist("[omega]"))
	{
		//get points without any transformation
		dr.getColumn(origValX, "[omega]");
		/*allocate space*/
		calcValX.resize(origValX.size());
		/*transform to radians*/
		transform(origValX.begin(), origValX.end(), calcValX.begin(),
		          bind2nd(std::multiplies<double>(), M_PI / 180.0));
		/*subtracts the x_max value from array*/
		toCenter(calcValX, origValZ);
	}
	else
	{
		throw Engine::Exception("Column \"[(d)omega]\" has not been found in "+
				m_programSettings->getEngineSettings().datafile);
	}

	/*allocate arguments and residuals*/
	for(size_t i = 0; i < calcValX.size(); ++i)
	{
		m_DataPoints.push_back(
				NonlinearFit::DataPoint(new CalculatorBaseArgument(calcValX[i]),
						origValZ[i]));
	}

	calcValZini.resize(origValZ.size(), 0.0);
	calcValZfin.resize(origValZ.size(), 0.0);
}

void Engine::toCenter(std::vector<double>& x, const std::vector<double>& y)
{
	std::vector<double>::const_iterator it_of_max;
	size_t index_of_max;
	double x_max/*, y_max*/;

	it_of_max = std::max_element(y.begin(), y.end());
	index_of_max = it_of_max - y.begin();

	x_max = x.at(index_of_max);
	//y_max = y.at(index_of_max);

	/*subtracts the x_max value from array*/
	transform(x.begin(), x.end(), x.begin(),
	          bind2nd(std::minus<double>(), x_max));

	//std::cout << "x_max:\t" << x_max << std::endl;
	//std::cout << "y_max:\t" << y_max << std::endl;
}

void Engine::setupCalculator()
{

	double Qx, Qz;
	const int * Q;
	double b_edge, b_screw;
	double lambda;
	double thickness;
	double nu;

	/*get Q in hexagonal miller indices*/
	Q = m_programSettings->getCalculatorSettings().Q;
	Qx = 2 * M_PI * sqrt(2.0 / 3 * (Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2]))
					  / m_programSettings->getSampleSettings().a0;
	Qz = 2 * M_PI * Q[3]
	                  / m_programSettings->getSampleSettings().c0;
	lambda = m_programSettings->getCalculatorSettings().lambda;

	b_edge = m_programSettings->getSampleSettings().b_edge;
	b_screw = m_programSettings->getSampleSettings().b_screw;

	thickness = m_programSettings->getSampleSettings().thickness;
	nu = m_programSettings->getSampleSettings().nu;

	try
	{
		switch (m_programSettings->getEngineSettings().workDiffractometry)
		{
		case ProgramSettings::EngineSettings::diffDOUBLE:
			m_calculator = new CalculatorDCRYST(Qx, Qz,
					b_edge, b_screw, lambda,
					nu, thickness);
			break;
		case ProgramSettings::EngineSettings::diffTRIPLE:
			m_calculator = new CalculatorTCRYST(Qx, Qz,
					b_edge, b_screw, lambda,
					nu, thickness);
			break;
		default:
			//this should never occur
			break;
		}
	} catch (std::exception& ex)
	{
		throw Engine::Exception(
				"Calculator allocation pb (" + toString(ex.what()) + ")");
	}
}

void Engine::addFParameter(const FitParameter & param)
{
	if ((m_programSettings->getEngineSettings().workMode
			== ProgramSettings::EngineSettings::modeFITI)
			&& (param.m_Lbvalue != param.m_Ubvalue))
	{
		m_fParametersInit.push_back(param);
	}
}

void Engine::setupCParameters()
{
	const FitParameter * param;
	/*scale*/
	param = &m_programSettings->getCalculatorSettings().scale;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*background*/
	param = &m_programSettings->getCalculatorSettings().background;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*densities*/
	/*edge*/
	param = &m_programSettings->getSampleSettings().edgeDislocations.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*screw*/
	param = &m_programSettings->getSampleSettings().screwDislocations.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*mixed*/
	param = &m_programSettings->getSampleSettings().mixedDislocations.rho;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*critical radii*/
	/*edge*/
	param = &m_programSettings->getSampleSettings().edgeDislocations.Rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*screw*/
	param = &m_programSettings->getSampleSettings().screwDislocations.Rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
	/*mixed*/
	param = &m_programSettings->getSampleSettings().mixedDislocations.Rc;
	m_cParameters[param->m_Name] = param->m_Value;
	addFParameter(*param);
}

void Engine::setupFitter()
{
	switch(m_programSettings->getEngineSettings().fitMethod)
	{
	case ProgramSettings::EngineSettings::fitLEVMAR:
		m_fitter = new LevmarFitter();
		break;
	/*case ProgramSettings::EngineSettings::fitNL2SOL:
	 * 	m_fitter = new Fitter(Fitter::ftNL2SOL, true);
		break;*/
	case ProgramSettings::EngineSettings::fitPORT:
	default:
		m_fitter = new PortFitter();
		break;
	}
	std::cout << "nb FitParameters:\t" << m_fParametersInit.size() << std::endl;

	m_fitter->init(m_calculator, m_fParametersInit, m_cParameters, m_DataPoints);
}

void Engine::printFitterInfo(NonlinearFit::NonlinearFitter * fitter)
{
	//size_t index;

	std::cout << "In " << fitter->getNbIterations() << " iterations ||f||^2 reduced from "
				<< fitter->getFinit() << " to " << fitter->getFfin() << std::endl;
	std::cout << "Nb Function evaluations:\t" << fitter->getNbFuncEval() << std::endl;
	std::cout << "Nb Jacobian evaluations:\t" << fitter->getNbGradEval() << std::endl;
	std::cout << "Reason to stop iterations:\t<" << fitter->getReasonToStop()
			<< "> Code:\t" << fitter->getReasonToStopId() << std::endl;

	std::cout << "Fitted parameters:" << std::endl;
	if(fitter->getName().compare("LEVMAR") == 0)
	{
		//index = 0;

		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << " +/- "
					<< sqrt(fitter->getCovarMatrix(i, i)) << std::endl;
		}

		std::cout << "Correlation matrix:" << std::endl;
		std::cout << "\t";
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
		}
		std::cout << std::endl;
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
			for(size_t j = 0; j <= i; ++j)
			{
				std::cout
						<< fitter->getCovarMatrix(i, j)
								/ sqrt(
										fitter->getCovarMatrix(i, i)
												* fitter->getCovarMatrix(j, j))
						<< "\t";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << std::endl;
		}
		std::cout << "Covariance matrix has not been calculated." << std::endl;
	}
}
