// lighttracer.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Photon.h"
#include "layer.h"
#include "material.h"
#include "Input.h"
#include "config.h"
#include "propagation.h"
#include "mcml_parser.h"
#include <iostream>
#include <fstream>
//Outdated-- will remove later...maybe
const float gran = 1000.0;
const double radius = 1;
//////////////////////////////////////
using namespace std;
struct dviwedi dvi;
struct skin ski;
struct renderoptions render;

int main()
{
	auto lay = layer(ski.refrac, ski.absorption, ski.scattering, ski.anistropy, propagation::getv0(ski.alpha));
	auto const material1 = material(render.photons, render.wth, 0.1, &lay, 0.999999999, false);

	//propagation prop(material1, res);
	
	//for (auto i = 0; i < prop.get_material().num_photons; i++)
	//{
	//	std::cout << 1- prop.get_material().matproperties->reflec << std::endl;
	//	photonstruct photon(1 - prop.get_material().matproperties->reflec, false, 0.0);
	//	prop.trace(photon); 
	//}
	//prop.write_to_logfile();
	

	auto mc = mcml::MCMLParser("C://Users//Tim//Documents//WISE17_18//Bachelor//FirstPrototype//mcml_test//bli.txt");
	std::cout << mc.GetEtaAbove() << std::endl;
	auto lay_info = mc.GetLayers();

	for(auto j : lay_info)
	{
		cout << j.mua_ << std::endl;
		cout << j.mus_ << "Scattering" << std::endl;
		cout << j.depth_ << "Tiefe" << std::endl;
		cout << j.eta_ << "Eta" << std::endl;
		cout << j.g_ << "anistropy" << std::endl;
	}
	//cin.get();
	return 0;
}




//Outdated-- will remove later...maybe
/*
 * This will update the buckets of the transmittance. It will check if there are already photons in a buck and
 * will add the new ones to it
 */
 //void updateBuckets(photonstruct &photon, std::vector<bucket> &buck)
 //{
 //	auto ex = false;
 //	auto const oncircle = line_circle(photon.position.x, photon.position.y);
 //	auto const x = round(oncircle.x * gran) / gran;
 //	auto const y = round(oncircle.y * gran) / gran;
 //	auto const z = round(photon.position.z * gran) / gran;
 //	for(auto i = 0; i < buck.size(); i++)
 //	{
 //		
 //		if(x == buck[i].x  && y == buck[i].y && z == buck[i].z)
 //		{
 //			buck[i].weight += photon.weight;
 //			buck[i].hits++;
 //			ex = true;
 //		}
 //		
 //	}
 //	if(!ex)
 //	{
 //		buck.push_back(bucket(x, y, z, photon.weight));
 //	}
 //}

 //inline double round(double const val)
 //{
 //	if (val < 0) return ceil(val - 0.5);
 //	return floor(val + 0.5);
 //}
 //
 //double specular_reflectance(layer &layer)
 //{
 //	//First value is set to one because the specular reflectance is 1 of vacuum
 //	return pow((1 - layer.refrac), 2) / pow(1 + layer.refrac, 2);
 //}