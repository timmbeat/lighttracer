#include "stdafx.h"
#include "dwivedi_sampling.h"
#include "config.h"
#include "constants.h"
#include "mcml_parser.h"

struct dviwedi dvi;
extern struct renderoptions render;




void dwivedi_sampling::run(const std::string mcml_path)
{
	//Parsing the Result of mcml
	auto mc = mcml::MCMLParser(mcml_path);
	auto const layer_mcml = mc.GetLayers();
	auto const layer_0 = layer_mcml[0];
 

	//Create layer and Material, these are the rendering options
	layer lay(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, getv0(layer_0.mus_ / (layer_0.mua_ + layer_0.mus_)));
	material material1(mc.get_numphotons(), render.wth, 0.1, &lay);
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
	write_to_logfile(&out, &material1, "dwivedi_logfile.txt", "dwivedi_output.csv");


}



double dwivedi_sampling::cal_stepsize(photonstruct* photon, material const * mat) 
{
	auto const y = random();
	auto const x = (mat->matproperties->v0 - 1) / (mat->matproperties->v0 + 1);
	auto const wz_new = mat->matproperties->v0*pow((mat->matproperties->v0 + 1), y)*(pow((mat->matproperties->v0 - 1), y) + pow((mat->matproperties->v0 + 1), y));
	auto const tmp = pow(x, y);
	auto const wz = mat->matproperties->v0 - (mat->matproperties->v0 + 1)*tmp;
	photon->wz = wz;
	return (-log(1 - random())) / ((1 - wz / mat->matproperties->v0)*(mat->matproperties->scattering + mat->matproperties->absorption));
}



double dwivedi_sampling::getv0(double const alpha)
{
	return alpha > 0.56 ? dvi.get_highv0(alpha) : dvi.get_lowv0(alpha);
}



void dwivedi_sampling::update_direction(photonstruct* photon, material const* mat)
{
	auto const wz = photon->wz;
	auto const y = 2 * slabProfiles::pi<double>() * random();
	auto const sint = sqrt(1 - wz * wz);
	auto const cosp = cos(y);
	double sinp;
	auto const ux = photon->direction.x;
	auto const uy = photon->direction.y;
	auto const uz = photon->direction.z;
	if (y < slabProfiles::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}

	if (fabs(photon->direction.z) > slabProfiles::cos_1<double>())
	{
		photon->direction.x = sint * cosp;
		photon->direction.y = sint * sinp;
		photon->direction.z = glm::sign(uz)*wz;
	}
	else
	{
		auto const temp = sqrt(1.0 - uz * uz);
		photon->direction.x = sint * (ux*uz*cosp - uy * sinp) / temp + ux * wz;

		photon->direction.y = sint * (uy*uz*cosp + ux * sinp) / temp + uy * wz;

		photon->direction.z = -sint * cosp*temp + uz * wz;
	}
}



