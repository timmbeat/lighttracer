#include "stdafx.h"
#include "dwivedi_sampling.h"
#include "config.h"
#include "constants.h"
#include "mcml_parser.h"
#include <iostream>
#include "gtc/constants.hpp"



#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/string_cast.hpp>
#include "gtx/norm.hpp"

#include "gtx/vector_angle.hpp"
#include "Logger.h"
#include <fstream>
#include "classic_sampling.h"
struct dviwedi dvi;
extern struct renderoptions render;
double dwivedi_distribution(double wz, double v0, double mu_t, double t);
double classical_distribution(double mu_t, double t);
double classical_stepsize(double t, double mu_t);
double anpassung(double wz, double v0);
double sampleThetaclas(double anisotropy);
void dwivedi_sampling::run(const std::string mcml_path)
{
	
	//Parsing the Result of mcml
	auto mc = mcml::MCMLParser(mcml_path);
	auto const layer_mcml = mc.GetLayers();
	auto const layer_0 = layer_mcml[0];
	classic_sampling clas{};
	Logger log{};
	//Create layer and Material, these are the rendering options
	layer lay(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, getv0(layer_0.mus_ / (layer_0.mua_ + layer_0.mus_)));
	const material material1(mc.get_numphotons(), render.wth, 0.1, &lay);
	
	//Create prop class for Rendering
	std::ofstream ccout("./plot/manyplots.csv", std::ofstream::trunc);
	std::stringstream csvout;
	csvout << std::setw(15) << std::left << "DWIVEDI" << std::setw(15) << std::left << "CLASSICAL" << "RUN";
	ccout << csvout.str() << std::endl;

	output out(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);
	output out_clas(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);

	csvout.str("");
	auto runs = 100;

	std::vector<double> sums;

	for (auto z = 0; z < runs; z++)
	{

		out = output(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);
		out_clas = output(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);

		for (auto i = 0; i < material1.num_photons; i++)
		{
			photonstruct photon(1 - lay.reflec, false);
			photonstruct photon_clas(1 - lay.reflec, false);
			while (!photon.dead)
			{
				trace(&photon, &out, &material1);
				if (photon.weight < material1.wth && !photon.dead)
				{

					roulette(&photon, &material1);
				}
		

			}
		
			while (!photon_clas.dead)
			{
				
				clas.trace(&photon_clas, &out_clas, &material1);
				if (photon_clas.weight < material1.wth && !photon_clas.dead)
				{

					roulette(&photon_clas, &material1);
				}
			
			}
	}

		auto dsum = 0.0;
		auto csum = 0.0;
		for (auto k = 0; k < out.bins_r; k++)
		{
			dsum += out.Rd_r[k][0];
		
			csum += out_clas.Rd_r[k][0];

		
		}
		//sums.push_back(dsum);
		csvout << std::setw(15) << std::left << dsum/material1.num_photons << std::setw(15) << std::left << csum/ material1.num_photons << z << " " << dsum/csum;
		ccout << csvout.str() <<std::endl;
		csvout.str("");
		std::cout << "RUN" << " " << z << std::endl;
	}
	//std::sort(sums.begin(), sums.end());

	//std::cout << abs(1.0-(sums[0] / sums[sums.size()-1]));
	log.create_PlotFile(mc, out, "plot/dwivedi_output.csv", material1);
	//log.create_RenderFile(mc, out_clas, "plot/logfile.txt", material1);
	log.create_RenderFile(mc, out, "plot/dwivedi_logfile.txt", material1);
	//write_to_logfile(&out, &material1, "dwivedi_logfile.txt", "dwivedi_output.csv", mc);


}



double dwivedi_sampling::cal_stepsize(photonstruct* photon, material const * mat) 
{
	
	auto const v0 = mat->matproperties->v0;
	auto const v0_1 = mat->matproperties->v0 + 1;

	//Immer negativ!
	auto const wz = v0 - v0_1 * pow((v0 - 1) / v0_1, random());

	photon->wz_old = photon->wz_new;
	photon->wz_new = wz;
	
	auto const ran = random();
	
	auto const dwistepsize = (-log(1 - ran)) /
		((1 - photon->wz_old /
		  mat->matproperties->v0)*
		  (mat->matproperties->scattering + mat->matproperties->absorption));

	
	///////////////////////////////////////////////////////////////////////////
	//////////////////////////Anisotropic Material/////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	auto const mu_t = mat->matproperties->mu_t;


	//This will yield, when isotropic materials given, a bad approximation to the true result. 
	//The curve will have strong Oszillation
	auto new_weight = photon->weight * classical_distribution(mu_t, dwistepsize) /
		dwivedi_distribution(photon->wz_old, v0, mu_t, dwistepsize);


	/*auto const new_weight = photon->weight * anpassung(wz, v0);
	*/
	//photon->weight = new_weight;

	////////////////////////////////////////////////////////////////////////////

	return dwistepsize;
}

double sampleThetaclas(double anisotropy)
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


	return scattering;
}


double anpassung(double wz, double v0)
{
	return (1 - wz / v0)*exp((1 - wz / v0));
}

double dwivedi_distribution(double wz, double v0, double mu_t, double t)
{
	
	//return exp(-(1 - wz_new / v0)*mu_t*t);

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
	auto wz = photon->wz_new;

	auto const g = mat->matproperties->anisotropy;
	auto const v0 = mat->matproperties->v0;
	//auto scattering = sampleThetaclas(g);

	//auto const theta = acos(-wz);
	auto const y = 2 * slabProfiles::pi<double>() * random();
	//auto const sint = sqrt(1 - wz * wz);
	auto const sint = sqrt(1-wz*wz);
    auto const cosp = cos(y);
	auto sinp = sqrt(1.0 - cosp*cosp);
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
			photon->direction.z = glm::sign(uz)*wz;


		auto const direction_new = photon->direction;
		auto const dot = glm::dot(direction_new, direction_old);
		auto const dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));
		photon->directiont = dot_theta;

		auto const cos_theta = dot_theta < 0.0 ? -1 - dot_theta : 1 - dot_theta;
		auto const hg = directional_distribution_hg(mat->matproperties->anisotropy, cos_theta); //Bis hier hin ist es richtig!
		auto const norm = glm::dvec3(0.0, 0.0, 1.0);
		auto const dot_norm = glm::dot(norm, direction_new);
		auto const dot_norm_theta = dot_norm / (glm::length(direction_new)*glm::length(norm));



		auto dwi = directional_distribution_dwi(mat->matproperties->v0, wz);
		


		
		auto new_weight = photon->weight * hg/dwi;

	
		photon->weight = new_weight;



}


double dwivedi_sampling::directional_distribution_dwi(double const v0, double const wz_)
{
	auto const res = (1 / log((v0 + 1) / (v0 - 1))) * (1 / (v0 - wz_));
	return res;

}

/*
 */
double dwivedi_sampling::directional_distribution_hg(double const g, double const theta)
{
	auto const qg = g * g;

	auto const res = (1 - qg) / (2 * pow(1 + qg - 2 * g * theta, 3.0 / 2.0));

	return res;
}

void dwivedi_sampling::cal_absorption(photonstruct* photon, material const* mat_) const
{

	auto path_mean = 1 / ((1 - photon->wz_old / mat_->matproperties->v0)*mat_->matproperties->mu_t);
	auto const tmp = photon->weight *mat_->matproperties->absorption  * path_mean;

	photon->weight -= tmp;
}

bool dwivedi_sampling::is_hit(photonstruct* photon, material const* mat_)
{

	if (photon->sleft == 0.0)
	{
		photon->step = cal_stepsize(photon, mat_);
	}
	else
	{
		photon->step = photon->sleft / ((1 - photon->wz_old / mat_->matproperties->v0)*mat_->matproperties->mu_t);
		photon->sleft = 0.0;
	}


	if (photon->direction.z < 0.0)
	{
		auto const s1 = -photon->position.z / photon->direction.z;
		//No check for lower boundery, because there is none
		if (s1 < photon->step && photon->direction.z != 0.0)
		{
			photon->sleft = (photon->step - s1)*((1 - photon->wz_old / mat_->matproperties->v0)*mat_->matproperties->mu_t);
			photon->step = s1;

			return true;
		}
	}

	return false;
}


