#include "stdafx.h"
#include "Logger.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "material.h"
#include <chrono>
#include <iostream>
#include <windows.h>
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

void Logger::create_RenderFolder(const mcml::MCMLParser& parser, output& values, const material& mat) const
{
	

	std::stringstream out;


	auto names = create_name_structure(parser.GetLayers()[0].mus_, parser.GetLayers()[0].mua_, parser.GetLayers()[0].g_, mat.num_photons);


	///////
	std::ofstream compout;
	compout.open(names[6], std::ofstream::out | std::ofstream::app);


	std::ofstream ccout(names[4], std::ofstream::trunc);
	auto mcml_ra = parser.get_ra();
	auto setw_nsize = get_longest_number(values.Rd_r) + 1;

	if (setw_nsize < SETW_SIZE) setw_nsize = SETW_SIZE;
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

			csvout << std::setw(setw_nsize) << std::left << values.Rd_r[j][i] * 1.0 / ((j + 0.5)*scale1) << std::setw(setw_nsize) << std::left << mcml_ra[j] << j << std::endl;
			ccout << csvout.str();
			counter += values.Rd_r[j][i];
		}
	}

	if (!values.sum_set)
	{
		values.sum = counter;
		values.sum_set = true;
	}
	std::stringstream compare_file;

	auto const diffuse_ref = counter / mat.num_photons;
	auto const diff = abs(parser.get_diffuse_reflectance() - diffuse_ref);
	compare_file << parser.get_diffuse_reflectance() << "," << diffuse_ref << "," << diff  << "," << abs(1-(parser.get_diffuse_reflectance()/(diffuse_ref))) * 100 << "," << mat.num_photons <<std::endl;
	compout << compare_file.str();
	compout.close();



	ccout.close();
	create_TeXFile(names[1], names[0]);


	std::stringstream systemcall;
	systemcall << "pdflatex " "-interaction=nonstopmode -output-directory=" << names[2] << " " << names[3];

	
	auto tmp = systemcall.str();

	const auto cstr = tmp.c_str();

	system(cstr);
	clean(names[1], names[0]);
}
/*
 * Delete the TeXFile and the mess what the compilation of the TeXfile created
 */
void Logger::clean(std::string file, std::string up_name) const
{
	std::stringstream loc;
	loc << ".\\mcml_examples\\Research\\" << up_name << "\\" << file << "\\" << file;

	std::stringstream syscall;

	syscall << "del " << loc.str() << ".log" << "," << loc.str() << ".aux" << "," << loc.str() << ".tex";
	auto tmp = syscall.str();
	syscall.str("");
	const auto cstr_2 = tmp.c_str();
	system(cstr_2);
}

/*
 * This Function will Create the TeXFile for creating the Plot
 */
void Logger::create_TeXFile(std::string name, std::string up_name) const
{
	
	std::stringstream locpath;
	std::stringstream csvPath;


	locpath << "./mcml_examples/Research/" << up_name << "/"<< name << "/" << name << ".tex";
	csvPath << "./mcml_examples/Research/" << up_name << "/"<< name << "/" << name << ".csv";

	std::ofstream texFile(locpath.str(), std::ofstream::trunc);
	
	std::stringstream tex;
	
	tex <<"\\documentclass{standalone} \n"
		  "\\usepackage[utf8]{inputenc} \n"
		  "\\usepackage{amsmath}\n"
		  "\\usepackage{amsfonts}\n"
		  "\\usepackage{amssymb}\n"
		  "\\usepackage{pgfplots}\n"
		  "\\usepackage{tikz,graphicx}\n"
		  "\\usetikzlibrary{arrows}\n"
		  "\\usetikzlibrary{external}\n"
		  "%\\tikzexternalize[prefix=./]\n"
		  "\\pgfplotsset{every tick label/.append style={font=\\tiny},\n"
		  "		every x label/.append style={yshift=.2em},\n"
		  "		every y label/.append style={xshift=.2em}}\n"
		  "\\newcommand{\\plotfilepath}{\"" << csvPath.str() << "\"}\n"
		  "\\begin{document}\n"
		  "\\begin{tikzpicture}\n"
		  "		\\begin{semilogyaxis}[width=0.95\\linewidth, height=5cm, xmin=0, xmax=100,\n"
		  "				  xlabel={$r$},\n"
		  "				  ylabel={$\\log\\left(R_d\\left(r\\right)\\right)$}]\n"
		  "			\\addplot[red] table[skip first n=0, mark=none,x index=2,y index=1]{\\plotfilepath};\n"
		  "			\\addplot[green] table[skip first n=0, mark=none, x index=2, y index=0]{\\plotfilepath};\n"
		  "		\\end{semilogyaxis}\n"
		  "\\end{tikzpicture}\n"
		  "\\end{document}";

	texFile << tex.str();
	texFile.close();
}

/*
 * Create a Folder 
 */
void Logger::CreateFolder(std::string name) const
{
	if(!CreateDirectory(s2ws(name).c_str(), NULL)) return;

	CreateDirectory(s2ws(name).c_str(), NULL);
}

std::string Logger::create_MCMLFile(double scattering, double absorption, double anisotropy, std::size_t num_photons)
{

	std::stringstream mcml;
	auto names = create_name_structure(scattering, absorption, anisotropy, num_photons);
	
	std::stringstream input_file;
	input_file << names[5] << "\\" << names[1] << "_input" <<".txt";
	std::ofstream mcmlfile(input_file.str(), std::ofstream::trunc);

	auto const in = input_file.str();


	input_file.str("");
	input_file << names[5] << "\\" << names[1] << "_result"<< ".txt";
	auto res = input_file.str();
	mcml << "1.0\n"
		"1\n"
		"" << input_file.str() << " A\n"
		<< num_photons << "\n" <<
		"0.0001 0.0001\n"
		"1 100 1\n"
		"1\n"
		"1.0\n"
		"1.3 " << absorption * 10 << " " << scattering * 10 << " " << anisotropy << " " << "1E+008\n"
		"1.0";

	mcmlfile << mcml.str() << std::endl;
	mcmlfile.close();

	std::stringstream syscall;

	syscall << "mcml " << in;

	auto tmp = syscall.str();
	auto const tmp_cstr = tmp.c_str();
	system(tmp_cstr);

	return res;
}

std::vector<std::string> Logger::create_name_structure(double scattering, double absorption, double anisotropy,
	std::size_t num_photons) const
{
	std::vector<std::string> names;


	std::stringstream up_name;
	std::stringstream name;
	std::stringstream loc;
	std::stringstream TeXFile;
	std::stringstream render_file;
	std::stringstream up_name_folder;
	std::stringstream compare_file;





	up_name << "abs" << absorption << "-scat" << scattering << "-g" << anisotropy;
	name << "pho" << num_photons << "-abs" << absorption << "-scat" << scattering << "-g" << anisotropy;
	loc << ".\\mcml_examples\\Research\\" << up_name.str() << "\\" << name.str();
	TeXFile << loc.str() << "\\" << name.str() << ".tex";
	render_file << loc.str() << "\\" << name.str() << ".csv";

	up_name_folder << ".\\mcml_examples\\Research\\" << up_name.str();
	compare_file << ".\\mcml_examples\\Research\\" << up_name.str() << "\\" << up_name.str() << ".txt";
	names.push_back(up_name.str());
	names.push_back(name.str());
	names.push_back(loc.str());
	names.push_back(TeXFile.str());
	names.push_back(render_file.str());
	names.push_back(up_name_folder.str());
	names.push_back(compare_file.str());
	CreateFolder(up_name_folder.str());

	CreateFolder(loc.str());

	return names;
}

/*
 * https://stackoverflow.com/questions/27220/how-to-convert-stdstring-to-lpcwstr-in-c-unicode
 */
std::wstring Logger::s2ws(const std::string& s) const
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}

/**
 * This Method will return the length of longest Number for formatting the CSV file
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
