#pragma once

#include <glm.hpp>

struct Layer
{
	Layer(double const refrac, double const absorption, double const scattering, double const anisotropy): 
	refrac(refrac), absorption(absorption), scattering(scattering), anisotropy(anisotropy)
	{
	}

	double refrac;
	double absorption;
	double scattering;
	double anisotropy;
};