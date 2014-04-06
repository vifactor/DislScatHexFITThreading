/*
 * main.cpp
 *
 *  Created on: 30 may 2012
 *      Author: kopp
 */

#include "Engine.h"
#include <iterator>

#define PROGRAM "DislScatThreadingVERT"
#define VERSION "1.55 from 4 april 2014"
#define AUTHOR "Viktor S. Kopp"

void info()
{
	std::cout <<"Program:\t" << PROGRAM << std::endl;
	std::cout <<"Version:\t" << VERSION << std::endl;
	std::cout <<"Author:\t" << AUTHOR << std::endl;
}

int main(int argc, char ** argv)
{
	std::string stgfile;
	ProgramSettings programSettings;
	Engine engine;

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

	std::cout << "Done successfully." << std::endl;
	return EXIT_SUCCESS;
}
