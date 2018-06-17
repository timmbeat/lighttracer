#pragma once

#include <glm.hpp>

struct layer
{
	layer(double const refrac, double const absorption, double const scattering, double const anisotropy, double const v0): 
	refrac(refrac), absorption(absorption), scattering(scattering), anisotropy(anisotropy),  v0(v0)
	{
		reflec = reflection(refrac);
	}

	double refrac;
	double absorption;
	double scattering;
	double anisotropy;
	double v0;
	double reflec;

private:
	static double reflection(double const refrac)
	{

		auto const tmp = (1 - refrac)/(1+refrac);
		return tmp * tmp;
	}
};