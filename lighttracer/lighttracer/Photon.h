#pragma once
#include <vector>
#include <glm.hpp>
#include <sstream>

struct Photonstruct
{
	glm::dvec3 position;
	glm::dvec3 directions;
	double weight;
	bool dead;
	//For now unnessary..
	short layer;
	double step;
	double stepleft;

	std::string to_string() const
	{
		std::stringstream output;
		output << "x " << position.x << " y " << position.y << " z " << position.z << std::endl;

		return output.str();
	}

};
