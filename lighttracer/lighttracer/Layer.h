#pragma once

/**
 * This structure will hold all information about the material which the photons moves in
 */
struct layer
{
	layer(double const refrac, double const absorption, double const scattering, double const anisotropy, double const v0): 
	refrac(refrac), absorption(absorption), scattering(scattering), anisotropy(anisotropy),  v0(v0)
	{
		reflec = reflection(refrac);
	}

	layer(double const refrac, double const absorption, double const scattering, double const anisotropy)
		:refrac(refrac), absorption(absorption), scattering(scattering), anisotropy(anisotropy)
	{
		reflec = reflection(refrac);
		v0 = 0;
	}
	double refrac;
	double absorption;
	double scattering;
	double anisotropy;
	double mu_t = scattering + absorption;
	double v0;
	double reflec;

private:
	static double reflection(double const refrac)
	{

		auto const tmp = (refrac - 1)/(refrac+1);
		
		return tmp * tmp;
	}
};