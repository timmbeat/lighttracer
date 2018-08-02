#pragma once
#include "mcml_parser.h"
#include "Input.h"

struct material;

class Logger
{
	public:
	Logger();
	~Logger();
	void create_PlotFile(const mcml::MCMLParser &parser, output &values, const std::string &plotfile, const material &mat) const;
	void create_RenderFile(const mcml::MCMLParser &parser, output &values, const std::string &plotfile, const material &mat) const;



	private:
	int get_longest_number(const std::vector < std::vector<double> > &numbers) const;
	std::string  getCurrentTimeAsString()const;
	int const SETW_SIZE = 15;
	int const FORMAT_30 = 30;
	int const FORMAT_5 = 5;

};

