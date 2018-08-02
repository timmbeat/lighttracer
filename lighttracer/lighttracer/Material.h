#pragma once
#include "layer.h"
#include <memory>

struct material
{
	material(long const num_photons, const double wth, const float rusroul, layer * matproperties): 
	num_photons(num_photons), wth(wth), rusroul(rusroul), matproperties(matproperties)
	{ }
	//the photons which will be traced
	const long num_photons;
	const double wth;
	const float rusroul;
	//the material in which we move
	layer *  matproperties;
	
};
