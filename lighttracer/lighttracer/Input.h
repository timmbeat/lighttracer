#pragma once
#include <glm.hpp>
#include <vector>
#include "bucket.h"
struct output
{
	output(int const bins_r, int const bins_a, double const delr, double const dela): bins_r(bins_r), bins_a(bins_a), delr(delr), dela(dela)
	{
		Rd_r = std::vector<std::vector<double>>(bins_r, std::vector<double>(bins_a));
	}
	std::vector<std::vector<double>> Rd_r;
	const int bins_r;
	const int bins_a;
	//double ref[ARRSIZEr][ARRSIZEa];
	double delr;
	double dela;
	std::vector<bucket> absorption;
	std::vector<bucket> trans;
};