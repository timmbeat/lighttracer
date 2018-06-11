#pragma once
#include <glm.hpp>
#include <vector>
#include "bucket.h"
struct output
{
	static const int ARRSIZE = 100;
	double ref[ARRSIZE][ARRSIZE];
	double delr = 0.1;
	double dela = 1;
	std::vector<bucket> absorption;
	std::vector<bucket> trans;
};