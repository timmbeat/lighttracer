#include "stdafx.h"
#include "Logger.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "material.h"
#include <chrono>
#include <iostream>

Logger::Logger() = default;


Logger::~Logger() = default;

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

void Logger::create_PlotFile(output& values_1, output& values_2, const std::string& plotfile, const material& mat) const
{

	//std::ofstream fout(path, std::ofstream::trunc);

	std::ofstream ccout(plotfile, std::ofstream::trunc);
	auto setw_nsize = get_longest_number(values_1.Rd_r) + 1;

	if (setw_nsize < SETW_SIZE) setw_nsize = SETW_SIZE;
	std::stringstream csvout;
	csvout << std::setw(setw_nsize) << std::left << "weight_1" << std::setw(setw_nsize) << std::left << "weight_2" << "ir";
	ccout << csvout.str() << std::endl;

	csvout.str("");
	std::stringstream output;



	auto counter = 0.0;
	auto counter_2 = 0.0;



	///////////////////
	//Scaling/////////
	//////////////////
	auto const scale1 = 2.0*slabProfiles::pi<double>()*values_1.delr*values_1.delr*mat.num_photons;


	for (auto j = 0; j < values_1.bins_r; j++)
	{
		for (auto i = 0; i < values_1.bins_a; i++)
		{
			output.str("");
			csvout.str("");

			csvout << std::setw(setw_nsize) << std::left << values_1.Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) << std::setw(setw_nsize) << std::left << values_2.Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) << j << std::endl;
			ccout << csvout.str();
			counter += values_1.Rd_r[j][i];
			counter_2 += values_2.Rd_r[j][i];
		}
	}

	if (!values_1.sum_set)
	{
		values_1.sum = counter;
		values_2.sum = counter_2;
		values_1.sum_set = true;
		values_2.sum_set = true;
	}
}

void Logger::create_RenderFile(const mcml::MCMLParser& parser, output& values, const std::string& renderfile,
							   const material& mat) const
{
	double transmittance;
	std::ofstream renderFile(renderfile, std::ofstream::trunc);
	std::stringstream out;



	if (!values.sum_set)
	{
		auto sum = 0.0;
		for (std::size_t i = 0; i < values.Rd_r.size(); i++)
		{
			for (std::size_t j = 0; j < values.Rd_r[0].size(); j++)
			{
				sum += values.Rd_r[i][j];
			}
		}
		values.sum = sum;
		values.sum_set = true;
		transmittance = sum / mat.num_photons;
	}

	else
	{
		transmittance = values.sum / mat.num_photons;
	}

	
	//Write the date to file 
	out << "IN THIS FILE YOU FIND INFORMATION ABOUT THE LAST SUCCESFUL RUN \n" <<
		"THIS RUN WAS SUCCESFUL ON " << getCurrentTimeAsString() << "\n\n" << std::endl;
	renderFile << out.str();
	out.str("");

	//Write to File the Transmittance of mcml and the Transmittance of the Run
	out << "MY TRANSMITTANCE " << std::setw(SETW_SIZE) << std::left << transmittance << "MCML TRANSMITTANCE " << std::left << parser.get_diffuse_reflectance() << "\n \n" <<std::endl;
	renderFile << out.str();
	out.str("");

	//Write all information about the material
	auto layers = parser.GetLayers();
	auto const layer = layers[0];

	out << std::setw(FORMAT_30) << std::left << "SCATTERING COEFFICIENT " << std::setw(FORMAT_5) << std::left << layer.mus_ << "\n"
		<< std::setw(FORMAT_30) << std::left << "ABSORPTION COEFFICIENT " << std::setw(FORMAT_5) << std::left << layer.mua_ << "\n"
		<< std::setw(FORMAT_30) << std::left << "ANISOTROPY " << std::setw(FORMAT_5) << std::left << layer.g_ << "\n" << std::endl;
	renderFile << out.str();
	out.str("");



}

void Logger::create_RenderFile(output& values_1, output& values_2, const std::string& renderfile,
	const material& mat) const
{
	double transmittance;
	double transmittance_2;
	std::ofstream renderFile(renderfile, std::ofstream::trunc);
	std::stringstream out;



	if (!values_1.sum_set || !values_2.sum_set)
	{
		auto sum = 0.0;
		auto sum_2 = 0.0;
		for (std::size_t i = 0; i < values_1.Rd_r.size(); i++)
		{
			for (std::size_t j = 0; j < values_1.Rd_r[0].size(); j++)
			{
				sum += values_1.Rd_r[i][j];
				sum_2 += values_2.Rd_r[i][j];
			}
		}
		values_1.sum = sum;
		values_2.sum = sum_2;
		values_1.sum_set = true;
		values_2.sum_set = true;
		transmittance = sum / mat.num_photons;
		transmittance_2 = sum_2 / mat.num_photons;
	}

	else
	{
		transmittance = values_1.sum / mat.num_photons;
		transmittance_2 = values_2.sum / mat.num_photons;
	}


	//Write the date to file 
	out << "IN THIS FILE YOU FIND INFORMATION ABOUT THE LAST SUCCESFUL RUN \n" <<
		"THIS RUN WAS SUCCESFUL ON " << getCurrentTimeAsString() << "\n\n" << std::endl;
	renderFile << out.str();
	out.str("");

	//Write to File the Transmittance of mcml and the Transmittance of the Run
	out << "VALUES_1 TRANSMITTANCE " << std::setw(SETW_SIZE) << std::left << transmittance << "VALUES_2 TRANSMITTANCE " << std::left << transmittance_2 << "\n \n" << std::endl;
	renderFile << out.str();
	out.str("");

	//Write all information about the material
	out << std::setw(FORMAT_30) << std::left << "SCATTERING COEFFICIENT " << std::setw(FORMAT_5) << std::left << mat.matproperties->scattering << "\n"
		<< std::setw(FORMAT_30) << std::left << "ABSORPTION COEFFICIENT " << std::setw(FORMAT_5) << std::left << mat.matproperties->absorption << "\n"
		<< std::setw(FORMAT_30) << std::left << "ANISOTROPY " << std::setw(FORMAT_5) << std::left << mat.matproperties->anisotropy << "\n" << std::endl;
	renderFile << out.str();
	out.str("");



}


/**
 * This Method will return the length of longest Number
 */
int Logger::get_longest_number(const std::vector<std::vector<double>>& numbers) const
{
	std::size_t max = 0;
	std::ostringstream strs;
	for (std::size_t i = 0; i < numbers.size(); i++)
	{
		for (std::size_t j = 0; j < numbers[0].size(); j++)
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
