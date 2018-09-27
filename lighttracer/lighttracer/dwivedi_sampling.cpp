/*
 * Created by Tim Mend for the Bachelor Thesis on Implementing the Dwivedi Sampling Scheme.
 */

#include "stdafx.h"
#include "dwivedi_sampling.h"
#include "config.h"
#include "constants.h"
#include "mcml_parser.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/string_cast.hpp>
#include "gtx/norm.hpp"
#include "Logger.h"
#include <thread>
#include "classic_sampling.h"
#include "output.h"
struct dviwedi dwi;

/*
 * This function provides a full 
 ing run for a given mcml example. 
 * It will first read a mcml result file with the MCMLParser provided by Sebastian Maisch
 * then creates the necassary structs for the rendering run(s). Then it will loops through
 * the amount of photons and will call the trace function. 
 */
void dwivedi_sampling::run(const std::string mcml_path)
{



	//Parsing the Result of mcml
	auto mc = mcml::MCMLParser(mcml_path);
	auto const layer_mcml = mc.GetLayers();
	auto const layer_0 = layer_mcml[0];
	Logger log{};

	//Create layer and Material, these are the rendering options
	layer lay(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, getv0(layer_0.mus_ / (layer_0.mua_ + layer_0.mus_)));
	const material material1(mc.get_numphotons(), get_threshould(), 0.1, &lay);
	output out(mc.get_numr(), mc.get_numa(), mc.get_dr_(), 1);
	for (auto i = 0; i < material1.num_photons; i++)
	{
		photonstruct photon(1 - lay.reflec, false);
		
			while(!photon.dead)
			{
				trace(&photon, &out, &material1);
				if(photon.weight < material1.wth && !photon.dead)
				{
					roulette(&photon, &material1);
				}
			}
	}


	log.create_PlotFile(mc, out, "mcml_examples/plot/dwivedi_output.csv", material1);
	log.create_RenderFile(mc, out, "mcml_examples/plot/dwivedi_logfile.txt", material1);

}
void dwivedi_sampling::run(double const mua, double const mus, double const anisotropy, std::size_t const photons, std::size_t const numr,
						   double const delr)
{
	Logger log{};
	classic_sampling clas{};
	layer lay(1.0, mua, mus, anisotropy, getv0(mus / (mua + mus)));
	const material material1(photons, get_threshould(), 0.1, &lay);
	output out(numr, 1, delr, 1);
	output out_clas(numr, 1, delr, 1);


	std::thread classical([&material1, &lay, &clas, &out_clas]()
	{
		for (auto i = 0; i < material1.num_photons; i++)
		{
			photonstruct photon_clas(1 - lay.reflec, false);

			while (!photon_clas.dead)
			{
				clas.trace(&photon_clas, &out_clas, &material1);
				if (photon_clas.weight < material1.wth && !photon_clas.dead)
				{
					clas.roulette(&photon_clas, &material1);
				}
			}
		}
	});

	std::thread dwivedi([&material1, &lay, this, &out]()
	{
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
	});

	classical.join();
	dwivedi.join();
	std::stringstream filename;
	std::stringstream filename_2;
	filename << "own_examples" << "/"<< "out-" << "abs" << mua << "-" << "scat" << mus << "-" << "anis" << anisotropy << "-" << "phot" << photons << ".csv";
	filename_2 << "own_examples" << "/" << "out-" << "abs" << mua << "-" << "scat" << mus << "-" << "anis" << anisotropy << "-" << "phot" << photons << ".txt";

	log.create_PlotFile(out_clas, out, filename.str(), material1);
	log.create_RenderFile(out_clas, out, filename_2.str(), material1);

}

/* 
 * This function will calculate the stepsize for the next move in the tissue. 
 * It will also normalize the weight for unbiased sampling.
 */
double dwivedi_sampling::cal_stepsize(photonstruct* photon, material const * mat)
{
	auto dwistepsize = 0.0;
	if(random() > mat->cpropability)
	{
		dwistepsize = sample_classical_distribution(mat->matproperties->mu_t);
	}

	else
	{
		dwistepsize = (-log(1 - random())) /
			((1 - photon->wz_old /
			  mat->matproperties->v0)*
			  (mat->matproperties->mu_t));


		auto const mu_t = mat->matproperties->mu_t;

		photon->weight = photon->weight * classical_path_distribution(mu_t, dwistepsize) /
			dwivedi_path_distribution(photon->wz_old, mat->matproperties->v0, mu_t, dwistepsize);
	}

	return dwistepsize;
}



/*
 * This function will get the right eigenvalue for dwivedi sampling
 */
double dwivedi_sampling::getv0(double const   alpha)
{
	return alpha > 0.56 ? dwi.get_highv0(alpha) : dwi.get_lowv0(alpha);
}



/*
 * This function will update the direction after a scattering event. 
 * It will also normalize the weight.
 */
void dwivedi_sampling::update_direction(photonstruct* photon, material const* mat)
{
	if (random() > mat->cpropability)
	{

		auto const v0 = mat->matproperties->v0;
		auto const v0_1 = mat->matproperties->v0 + 1;

		auto const wz = v0 - v0_1 * pow((v0 - 1) / v0_1, random());

		photon->wz_old = photon->wz_new;
		photon->wz_new = wz;

		auto const y = 2 * slabProfiles::pi<double>() * random();
		auto const sint = sqrt(1 - wz * wz);
		auto const cosp = cos(y);
		double sinp;
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

		photon->direction = glm::normalize(photon->direction);


		auto const direction_new = photon->direction;
		auto const cos_theta = glm::dot(direction_new, direction_old);




		auto const hg = directional_distribution_hg(mat->matproperties->anisotropy, acos(cos_theta));
		auto const dwi = directional_distribution_dwi(mat->matproperties->v0, wz);


		photon->weight = photon->weight * hg / dwi;
	}

	else
	{
		propagation::update_direction(photon, mat);
	}



}

/*
 * This is the dwivedi pdf for the directional distribution
 */
double dwivedi_sampling::directional_distribution_dwi(double const v0, double const wz) const
{
	return (1 / log((v0 + 1) / (v0 - 1))) * (1 / (v0 - wz));
}

/*
 * This is the henyey Greenstein function 
 */
double dwivedi_sampling::directional_distribution_hg(double const g, double const theta) const
{
	auto const qg = g * g;

	return (1.0 - qg) / (2 * pow(1.0 + qg - 2.0 * g * cos(theta), 3.0 / 2.0));
}

/*
 * Calculates the absorption
 */
void dwivedi_sampling::cal_absorption(photonstruct* photon, material const* mat) const
{
	auto const tmp = photon->weight * mat->matproperties->absorption / (mat->matproperties->mu_t);
	/*auto path_mean = 1 / ((1 - photon->wz_old / mat->matproperties->v0)*mat->matproperties->mu_t);
	auto const tmp = photon->weight *mat->matproperties->absorption  * path_mean;
	*/
	photon->weight -= tmp;
}

/*
 * This is the pdf for the path distribution
 */
double dwivedi_sampling::dwivedi_path_distribution(double const wz, double const v0, double const mu_t, double const t) const
{
	return (1 - wz / v0)*mu_t*exp(-(1 - wz / v0)*mu_t*t);
}



