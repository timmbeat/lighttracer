#include "stdafx.h"
#include "classic_sampling.h"
#include <random>
#include "constants.h"
#include "mcml_parser.h"
#include "config.h"

extern struct renderoptions render;
classic_sampling::classic_sampling()
= default;


classic_sampling::~classic_sampling()
= default;


void classic_sampling::run(const std::string mcml_path)
{
	//Parsing the Result of mcml
	auto mc = mcml::MCMLParser(mcml_path);
	auto const layer_mcml = mc.GetLayers();
	auto const layer_0 = layer_mcml[0];


	//Create layer and Material, these are the rendering options
	layer lay(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, 0);
	material material1(mc.get_numphotons(), render.wth, 0.1, &lay, 0.999999999, false);
	output out(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);


	//Create prop class for Rendering
	for (auto i = 0; i < material1.num_photons; i++)
	{
		photonstruct photon(1 - lay.reflec, false);
		while (!photon.dead)
		{
			trace(&photon, &out, &material1);
			if (photon.weight < material1.wth && !photon.dead)
			{
				roulette(&photon, &material1);
			}
		}
	}
	write_to_logfile(&out, &material1);


}



double classic_sampling::cal_stepsize(photonstruct* photon, material const* mat)
{
	double ran;
	do
	{
		ran = random();

	}
	while (ran <= 0.0);
	auto const ret = -log(ran) / (mat->matproperties->absorption + mat->matproperties->scattering);

	return ret;
}




glm::dvec2 classic_sampling::calculate_scattering(double anisotropy)
{
	double scattering;
	if (anisotropy == 0.0)
	{

		scattering = 2 * random() - 1;
	}
	else
	{
		auto const g = anisotropy;
		auto const temp = (1 - g * g) / (1 - g + 2 * g * random());
		scattering = (1 + g * g - temp * temp) / (2 * g);
	}

	auto const azimuthal = 2 * slabProfiles::pi<double>() * random();

	return glm::dvec2(scattering, azimuthal);
}


void classic_sampling::update_direction(photonstruct* photon, material const* mat)
{
	auto const angles = calculate_scattering(mat->matproperties->anisotropy);
	//auto const angles = calculate_scattering();
	auto const sint = sqrt(1 - angles.x * angles.x);
	auto const cosp = cos(angles.y);
	double sinp;
	auto const ux = photon->direction.x;
	auto const uy = photon->direction.y;
	auto const uz = photon->direction.z;
	if (angles.y < slabProfiles::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}

	if (fabs(photon->direction.z) > 0.99999999)
	{
		photon->direction.x = sint * cosp;
		photon->direction.y = sint * sinp;
		photon->direction.z = glm::sign(uz)*angles.x;
	}
	else
	{
		auto const temp = sqrt(1.0 - uz * uz);
		photon->direction.x = sint * (ux*uz*cosp - uy * sinp) / temp + ux * angles.x;

		photon->direction.y = sint * (uy*uz*cosp + ux * sinp) / temp + uy * angles.x;

		photon->direction.z = -sint * cosp*temp + uz * angles.x;
	}
}

