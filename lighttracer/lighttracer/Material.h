#pragma once
#include "layer.h"

struct material
{
	material(long const num_photons, const double wth, const float rusroul, layer * matproperties, const double angle_threshold, const bool dwivedi): 
	num_photons(num_photons), wth(wth), rusroul(rusroul), angle_threshold(angle_threshold), dwivedi(dwivedi), matproperties(matproperties)
	{ }
	//the photons which will be traced
	const long num_photons;
	const double wth;
	const float rusroul;
	const double angle_threshold;
	const bool dwivedi;
	//the material in which we move
	layer *matproperties;

	
};