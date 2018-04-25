#pragma once
#include "Layer.h"

struct material
{
	material(long const num_photons, const double wth, Layer * matproperties): 
	num_photons(num_photons), wth(wth), matproperties(matproperties)
	{ }
	//the photons which will be traced
	const long num_photons;
	const double wth;
	//the material in which we move
	Layer *matproperties;
};