#include "stdafx.h"
#include "Logger.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "material.h"
#include <chrono>
#include <iostream>

Logger::Logger()
{
}


Logger::~Logger()
{
}

void Logger::create_PlotFile(const mcml::MCMLParser& parser, output& values, const std::string& plotfile, const material &mat) const
{

	//std::ofstream fout(path, std::ofstream::trunc);

	std::ofstream ccout(plotfile, std::ofstream::trunc);
	auto mcml_ra = parser.get_ra();
	auto setw_nsize = get_longest_number(values.Rd_r) + 1;

	if (setw_nsize < SETW_SIZE) setw_nsize = SETW_SIZE ;
	std::stringstream csvout;
	csvout << std::setw(setw_nsize) << std::left << "weight_me" << std::setw(setw_nsize) << std::left << "weight_mcml" << "ir";
	ccout << csvout.str() << std::endl;

	csvout.str("");
	std::stringstream output;



	auto counter = 0.0;




	///////////////////
	//Scaling/////////
	//////////////////
	auto const scale1 = 2.0*slabProfiles::pi<double>()*values.delr*values.delr*mat.num_photons;


	for (auto j = 0; j < values.bins_r; j++)
	{
		for (auto i = 0; i < values.bins_a; i++)
		{
			output.str("");
			csvout.str("");

			csvout << std::setw(setw_nsize) << std::left << values.Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) <<  std::setw(setw_nsize) << std::left << mcml_ra[j] << j << std::endl;
			ccout << csvout.str();
			counter += values.Rd_r[j][i];
		}
	}

	if (!values.sum_set)
	{
		values.sum = counter;
		values.sum_set = true;
	}

}

void Logger::create_RenderFile(const mcml::MCMLParser& parser, output& values, const std::string& renderfile,
							   const material& mat) const
{
	auto Transmittance = 0.0;
	std::ofstream renderFile(renderfile, std::ofstream::trunc);
	std::stringstream out;



	if (!values.sum_set)
	{
		auto sum = 0.0;
		for (auto i = 0; i < values.Rd_r.size(); i++)
		{
			for (auto j = 0; j < values.Rd_r[0].size(); j++)
			{
				sum += values.Rd_r[i][j];
			}
		}
		values.sum = sum;
		values.sum_set = true;
		Transmittance = sum / mat.num_photons;
	}

	else
	{
		Transmittance = values.sum / mat.num_photons;
	}

	
	//Write the date to file 
	out << "IN THIS FILE YOU FIND INFORMATION ABOUT THE LAST SUCCESFUL RUN \n" <<
		"THIS RUN WAS SUCCESFUL ON " << getCurrentTimeAsString() << "\n\n" << std::endl;
	renderFile << out.str();
	out.str("");

	//Write to File the Transmittance of mcml and the Transmittance of the Run
	out << "MY TRANSMITTANCE " << std::setw(SETW_SIZE) << std::left << Transmittance << "MCML TRANSMITTANCE " << std::left << parser.get_diffuse_reflectance() << "\n \n" <<std::endl;
	renderFile << out.str();
	out.str("");

	//Write all information about the material
	auto layers = parser.GetLayers();
	auto layer = layers[0];

	out << std::setw(FORMAT_30) << std::left << "SCATTERING COEFFICIENT " << std::setw(FORMAT_5) << std::left << layer.mus_ << "\n"
		<< std::setw(FORMAT_30) << std::left << "ABSORPTION COEFFICIENT " << std::setw(FORMAT_5) << std::left << layer.mua_ << "\n"
		<< std::setw(FORMAT_30) << std::left << "ANISOTROPY " << std::setw(FORMAT_5) << std::left << layer.g_ << "\n" << std::endl;
	renderFile << out.str();
	out.str("");



}


/**
 * This Method will return the length of longest Number
 */
int Logger::get_longest_number(const std::vector<std::vector<double>>& numbers) const
{
	auto max = 0;
	std::ostringstream strs;
	for (auto i = 0; i < numbers.size(); i++)
	{
		for (auto j = 0; j < numbers[0].size(); j++)
		{
			strs.str("");
			strs << numbers[i][j];
			auto number_str = strs.str();
			auto const length = number_str.size();
			if (length > max) max = length;

		}
	}
	return max;
}


/*
 * This returns the current Time without Milliseconds
 */
std::string Logger::getCurrentTimeAsString() const
{
	char str[26];
	auto const time_t = std::chrono::system_clock::now();
	const auto t = std::chrono::system_clock::to_time_t(time_t);
	ctime_s(str, sizeof str, &t);

	return str;
}
