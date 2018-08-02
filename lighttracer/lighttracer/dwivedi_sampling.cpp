#include "stdafx.h"
#include "dwivedi_sampling.h"
#include "config.h"
#include "constants.h"
#include "mcml_parser.h"
#include <iostream>
#include "gtc/constants.hpp"

#define GLM_ENABLE_EXPERIMENTAL

#include "gtx/norm.hpp"

#include "gtx/vector_angle.hpp"
#include "Logger.h"
struct dviwedi dvi;
extern struct renderoptions render;
double dwivedi_distribution(double wz, double v0, double mu_t, double t);
double classical_distribution(double mu_t, double t);
double classical_stepsize(double t, double mu_t);
void dwivedi_sampling::run(const std::string mcml_path)
{
	//Parsing the Result of mcml
	auto mc = mcml::MCMLParser(mcml_path);
	auto const layer_mcml = mc.GetLayers();
	auto const layer_0 = layer_mcml[0];
 
	Logger log{};
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

	log.create_PlotFile(mc, out, "dwivedi_output.csv", material1);
	log.create_RenderFile(mc, out, "dwivedi_logfile.txt", material1);
	//write_to_logfile(&out, &material1, "dwivedi_logfile.txt", "dwivedi_output.csv", mc);


}



double dwivedi_sampling::cal_stepsize(photonstruct* photon, material const * mat) 
{
	
	auto const v0 = mat->matproperties->v0;
	auto const v0_1 = mat->matproperties->v0 + 1;

	//Immer negativ!
	auto const wz = v0 - v0_1 * pow((v0 - 1) / v0_1, random());
	
	photon->wz = wz;
	
	auto const ran = random();
	
	auto const dwistepsize = (-log(1 - ran)) /
		((1 - wz /
		  mat->matproperties->v0)*
		  (mat->matproperties->scattering + mat->matproperties->absorption));

	
	///////////////////////////////////////////////////////////////////////////
	//////////////////////////Anisotropic Material/////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	auto const mu_t = mat->matproperties->mu_t;


	//This will yield, when isotropic materials given, a bad approximation to the true result. 
	//The curve will have strong Oszillation
	auto const new_weight = photon->weight * classical_distribution(mu_t, dwistepsize) /
		dwivedi_distribution(wz, v0, mu_t, dwistepsize);
	   

	
	photon->weight = new_weight;

	////////////////////////////////////////////////////////////////////////////

	return dwistepsize;
}

double dwivedi_distribution(double wz, double v0, double mu_t, double t)
{
	
	//return exp(-(1 - wz / v0)*mu_t*t);

	return (1 - wz / v0)*mu_t*exp(-(1 - wz / v0)*mu_t*t);

}

double classical_distribution(double const mu_t, double const t)
{
	return mu_t*exp(-mu_t * t);
}

double dwivedi_sampling::getv0(double const   alpha)
{
	return alpha > 0.56 ? dvi.get_highv0(alpha) : dvi.get_lowv0(alpha);
}


double classical_stepsize(double const t, double const mu_t)
{
	return -log(t) / mu_t;
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
	auto const direction_old = photon->direction;
	if (y < slabProfiles::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}

		photon->direction.x = sint * cosp;
		photon->direction.y = sint * sinp;
		photon->direction.z = glm::sign(uz)* wz;


	
		auto const direction_new = photon->direction;
		auto const dot = glm::dot(direction_new, direction_old);
		auto const dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));
		auto const hg = scattering_function_hg(mat->matproperties->anisotropy, dot_theta); //Bis hier hin ist es richtig!
	

		auto dwi = scattering_function(mat->matproperties->v0, wz);



		auto tmp = hg / dwi;
	
		auto new_weight = photon->weight * tmp;


		photon->weight = new_weight;


		//if (photon->weight > 1 || photon->weight == NAN)
		//{
		//	if (photon->weight == NAN) std::cout << "ERROR::PHOTON WEIGHT = NAN" << std::endl;

		//	else
		//	{
		//		std::cout << "ERROR:::PHOTON TO HEAVY " << photon->weight << std::endl;
		//		std::cout << "HENVEY GREENSTEIN::: " << hg << "  DWIVEDI::: " << dwi <<std::endl;
		//		std::cout << tmp << std::endl;
		//		std::cin.get();
		//	}
		//}
		//if(photon->weight > 1.0 || hg < 0.0)
		//{
		//	std::cout << "HG " << hg << " DWI " << dwi << std::endl;
		//	//std::cin.get();
		//}

		//photon->weight = photon->weight *  scattering_function_hg(mat->matproperties->anisotropy, uz) / wz;


}


double dwivedi_sampling::scattering_function(double const v0, double const wz_)
{
	auto const res = (1 / log((v0 + 1) / (v0 - 1))) * (1 / (v0 - wz_));
	return res;

}

double dwivedi_sampling::scattering_function_hg(double const g, double const theta)
{
	auto const qg = g * g;

	auto const res = (1 - qg) / (2 * pow(1 + qg - 2 * g * theta, 3 / 2));

	return res;
}

void dwivedi_sampling::cal_absorption(photonstruct* photon, material const* mat_) const
{
	auto const tmp = photon->weight * mat_->matproperties->absorption / ((1-photon->wz/mat_->matproperties->v0)*mat_->matproperties->mu_t);
	photon->weight -= tmp;
}

bool dwivedi_sampling::is_hit(photonstruct* photon, material const* mat_)
{
	//TODO: Test for both sampling sampling_scheme

	if (photon->sleft == 0.0)
	{
		photon->step = cal_stepsize(photon, mat_);
	}
	else
	{
		photon->step = photon->sleft / ((1 - photon->wz / mat_->matproperties->v0)*mat_->matproperties->mu_t);
		photon->sleft = 0.0;
	}


	if (photon->direction.z < 0.0)
	{
		auto const s1 = -photon->position.z / photon->direction.z;
		//No check for lower boundery, because there is none
		if (s1 < photon->step && photon->direction.z != 0.0)
		{
			photon->sleft = (photon->step - s1)*((1 - photon->wz / mat_->matproperties->v0)*mat_->matproperties->mu_t);
			photon->step = s1;

			return true;
		}
	}

	return false;
}


