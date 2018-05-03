#pragma once
#include <glm.hpp>
#include <vector>
struct InOut
{
	
	std::vector<glm::vec3> tissue;
	std::vector<glm::vec3> absorption;
	std::vector<Photonstruct> transmittance;

};