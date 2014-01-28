/*
 * main.cpp
 *
 *  Created on: 30 may 2012
 *      Author: kopp
 */

#include "Engine.h"
#include <Log.h>
#include <iterator>

#define PROGRAM "DislScatThreadingVERT"
#define VERSION "1.54 from 10 october 2013"
#define AUTHOR "Viktor S. Kopp"

void info()
{
	LOG(logTIME) << std::endl;
	LOG(logINFO)<<"Program:\t"<<PROGRAM << std::endl;
	LOG(logINFO)<<"Version:\t"<<VERSION << std::endl;
	LOG(logINFO)<<"Author:\t"<<AUTHOR << std::endl;
}

int main(int argc, char ** argv)
{
	std::string stgfile;
	ProgramSettings programSettings;
	Engine engine;

	OutputTO::Stream()=stdout;
	info();
	if(argc == 1)
	{
		stgfile = "default.cfg";
		try
		{
			programSettings.read(stgfile);
			engine.exec(&programSettings);
		}catch(ProgramSettings::Exception& ex)
		{
			std::cout << ex.what() << std::endl;
		}catch(Engine::Exception& ex)
		{
			std::cout << ex.what() << std::endl;
		}

	}
	else
	{
		for (int iarg = 1; iarg < argc; ++iarg)
		{
			stgfile = argv[iarg];
			try
			{
				programSettings.read(stgfile);
				engine.exec(&programSettings);
			}catch(ProgramSettings::Exception& ex)
			{
				std::cout << ex.what() << std::endl;
			}catch(Engine::Exception& ex)
			{
				std::cout << ex.what() << std::endl;
			}
		}
	}

	OutputTO::Stream()=stdout;
	std::cout << "Done successfully." << std::endl;
	return EXIT_SUCCESS;
}
