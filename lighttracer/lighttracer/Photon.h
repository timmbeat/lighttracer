#pragma once
#include <vector>
#include <glm.hpp>
#include <sstream>

struct photonstruct
{
	photonstruct(double const weight, bool const dead) :
		position(glm::dvec3(0, 0, 0)), direction(glm::dvec3(0, 0, 1)), weight(weight), dead(dead), sleft(0), alreadymoved(false)
	{
	}
	glm::dvec3 position;
	glm::dvec3 direction;
	double weight;
	int moves = 0;
	bool dead;
	//experimental
	bool alreadymoved;
	//For now unnessary..
	double step = 0.0;
	double sleft = 0.0;
	std::string to_string() const
	{
		std::stringstream output;
		output << "x " << direction.x << " y " << direction.y << " z " << direction.z << std::endl;

		return output.str();
	}

};
