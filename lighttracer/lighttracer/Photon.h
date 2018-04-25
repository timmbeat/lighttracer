#pragma once
#include <vector>
#include <glm.hpp>

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

};