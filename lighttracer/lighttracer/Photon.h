#pragma once
#include <vector>
#include "glm/glm.hpp"
#include <sstream>
#include <iomanip>

struct photonstruct
{
	photonstruct(double const weight, bool const dead) :
		position(glm::dvec3(0, 0, 0)), direction(glm::dvec3(0, 0, 1)), weight(weight), dead(dead), sleft(0)
	{
	}
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;
	bool dead;
	//For now unnessary..
	double step = 0.0;
	double sleft = 0.0;
	double wz_old = -1.0;
	std::string to_string() const
	{
		std::stringstream output;
		output  << std::setw(15) << std::left << "DIRECTION" << std::setw(15) << std::left << direction.x << std::setw(15) << std::left << direction.y << std::setw(15) << std::left <<  direction.z << "\n"
				<< std::setw(15) << std::left << "POSITION"  << std::setw(15) << std::left << position.x  << std::setw(15) << std::left << position.y  << std::setw(15) << std::left <<  position.z;

		return output.str();
	}

};
