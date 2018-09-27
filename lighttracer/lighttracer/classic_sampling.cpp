#include "stdafx.h"
#include "classic_sampling.h"
#include "constants.h"
#include "mcml_parser.h"
#include "Logger.h"

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

	Logger log{};

	//Create layer and Material, these are the rendering options
	layer lay(layer_0.eta_, layer_0.mua_, layer_0.mus_, layer_0.g_, 0);
	material material1(mc.get_numphotons(), get_threshould(), 0.1, &lay);
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
	log.create_RenderFile(mc, out, "logfile.txt", material1);
	log.create_PlotFile(mc, out, "output.csv", material1);

}

double classic_sampling::cal_stepsize(photonstruct* photon, material const* mat)
{
	return sample_classical_distribution(mat->matproperties->mu_t);
}
