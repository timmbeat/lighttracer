#pragma once
#include "layer.h"
#include <memory>
#include <algorithm>

struct material
{
	material(long const num_photons, const double wth, const double rusroul, layer * matproperties): 
	num_photons(num_photons), wth(wth), rusroul(rusroul), matproperties(matproperties), cpropability(std::max(0.1, pow(abs(matproperties->anisotropy), 1.0 / 3.0)))
	{ }
	//the photons which will be traced
	const long num_photons;
	const double wth;
	const double rusroul;
	const double cpropability;
	//the material in which we move
	layer *  matproperties;
	
};
