#pragma once
#include "mcml_parser.h"
#include "output.h"

struct material;

class Logger
{
	public:
	Logger();
	~Logger();
	void create_PlotFile(const mcml::MCMLParser &parser, output &values, const std::string &plotfile, const material &mat) const;
	void create_PlotFile(output &values_1, output &values_2, const std::string &plotfile, const material &mat) const;
	void create_RenderFile(const mcml::MCMLParser &parser, output &values, const std::string &plotfile, const material &mat) const;
	void create_RenderFile(output &values_1, output &values_2, const std::string &plotfile, const material &mat) const;
	void create_RenderFolder(const mcml::MCMLParser &parser, output &values, const material &mat)const;
	void create_TeXFile(std::string name, std::string up_name) const;
	void CreateFolder(std::string name) const;
	std::string create_MCMLFile(double scattering, double absorption, double anisotropy, std::size_t num_photons);

	private:
	std::vector<std::string> create_name_structure(double scattering, double absorption, double anisotropy, std::size_t num_photons) const;
	int get_longest_number(const std::vector < std::vector<double> > &numbers) const;
	std::string  getCurrentTimeAsString()const;
	int const SETW_SIZE = 15;
	int const FORMAT_30 = 30;
	int const FORMAT_5 = 5;



	private:
	std::wstring s2ws(const std::string& s) const;
	void clean(std::string file, std::string up_name) const ;
};

