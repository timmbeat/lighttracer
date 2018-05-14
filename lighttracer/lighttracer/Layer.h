#pragma once

#include <glm.hpp>

struct Layer
{
	Layer(double const refrac, double const absorption, double const scattering, double const anisotropy, double const v0): 
	refrac(refrac), absorption(absorption), scattering(scattering), anisotropy(anisotropy),  v0(v0)
	{
	}

	double refrac;
	double absorption;
	double scattering;
	double anisotropy;
	double v0;
};