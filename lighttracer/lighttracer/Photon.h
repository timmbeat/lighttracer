#pragma once
#include <vector>
#include <glm.hpp>
#include <sstream>

struct photonstruct
{
	photonstruct(double const weight, bool const dead, double const step) :
		position(glm::dvec3(0, 0, 0)), direction(glm::dvec3(0, 0, 1)), weight(weight), dead(dead), step(step), sleft(0)
	{
	}
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;
	bool dead;
	//For now unnessary..
	double step;
	double sleft;
	std::string to_string() const
	{
		std::stringstream output;
		output << "x " << position.x << " y " << position.y << " z " << position.z << std::endl;

		return output.str();
	}

};
