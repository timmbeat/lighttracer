// lighttracer.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Photon.h"
#include "Layer.h"
#include "material.h"
#include "Input.h"
#include <random>
#include <iostream>
#include <math.h>
#include <vector>
#include <sstream>
#include <string>
#include "config.h"
#include <complex>
#include <fstream>
#include "bucket.h";

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);
const double angle_threshold = 0.99999;
const double m_pi = 3.14159265358979323846;
const double chance = 0.1;
const float gran = 1000.0;
const double radius = 1;
InOut res;
std::string sampling_scheme = "";

///////////////////////////////////////////////////////////////////////////////////////////
//Set this variable to true if you want to use dwivedi sampling
///////////////////////////////////////////////////////////////////////////////////////////
bool dwivedi_sam = false;

void absorption(material &material, Photonstruct &photon);
double step(Layer &layer);
glm::dvec2 calculate_scattering(material &material);
void move(Photonstruct &photon, double stepsize);
void update_direction(Photonstruct &photon, material &material);
double specular_reflectance(Layer &layer);
double sign(double x);
void trace(material &mat, Photonstruct &photon);
void render(material &mat);
double getv0(double const alpha);
void updateBuckets(Photonstruct photon);
glm::dvec2 line_circle(double x, double y);

void write_text_to_log_file(const std::string &text);

int main()
{
	if(dwivedi_sam)
	{
		sampling_scheme = " DWIVDI ";
	}
	else
	{
		sampling_scheme = " STANDART ";
	}
	for(auto iter = 0; iter < 1; iter++)
	{
		
		std::cout << iter << std::endl;

	auto layer = Layer(skin.absorption, skin.scattering, skin.anistropy, skin.refrac, getv0(skin.alpha));
	auto material1 = material(renderoptions.photons, renderoptions.wth, &layer);
	render(material1);

	auto counter = 0.0;
	auto poscountx = 0.0;
	auto poscounty = 0.0;
	for (const auto i : res.transmittance)
	{
		counter += i.weight;
		poscountx += abs(i.position.x);
		poscounty += abs(i.position.y);
		
	}
	std::stringstream buckets;
	std::sort(res.trans.begin(), res.trans.end(), sort_bucket());
	for(const auto p : res.trans)
	{
		std::cout << p.weight << std::endl;
		buckets.str("");
		buckets << "WEIGHT " << p.weight << " POSITION X " << p.x << " Y " << p.y << " Z " << p.z << " HITS " << p.hits;
		write_text_to_log_file(buckets.str());
	}
	std::stringstream output;
	
	//output << "WEIGHT " << counter << " DISTANCE X " << poscountx << " DISTANCE Y " << poscounty << " NORMAL SAMPLING";

	output << "WEIGHT " << counter/renderoptions.photons << " DISTANCE X " << poscountx << " DISTANCE Y " << poscounty << sampling_scheme;
	
		write_text_to_log_file(output.str());
	/*std::cout << "Absorption" << std::endl;
	for (const auto i : res.absorption)
	{
		
		std::cout << "x " << i.x << " y " << i.y << " z " << i.z << std::endl;
	}
	*/

	
	}
	return 0;
}

void absorption(material& material, Photonstruct& photon)
{
	auto const tmp = photon.weight * (material.matproperties->absorption/(material.matproperties->absorption + material.matproperties->scattering));
	photon.weight -= tmp;
	//push absorption into infinite array
	res.absorption.push_back(glm::vec3(photon.position.x, photon.position.y, photon.position.z));
}

void write_text_to_log_file(const std::string &text)
{
	std::ofstream fout("logfile.txt", std::ofstream::app);
	fout << text;
	fout << std::endl;
}

void updateBuckets(Photonstruct photon)
{
	auto ex = false;
	auto const oncircle = line_circle(photon.position.x, photon.position.y);
	auto const x = round(oncircle.x * gran) / gran;
	auto const y = round(oncircle.y * gran) / gran;
	auto const z = round(photon.position.z * gran) / gran;
	for(auto i = 0; i < res.trans.size(); i++)
	{
		
		if(x == res.trans[i].x  && y == res.trans[i].y && z == res.trans[i].z)
		{
			std::cout << "!!!!!!!!!!!!!!!HIT" << std::endl;
			res.trans[i].weight += photon.weight;
			res.trans[i].hits++;
			ex = true;
		}
		
	}
	if(!ex)
	{
		res.trans.push_back(bucket(x, y, z, photon.weight));
	}
	
}
inline double round(double val)
{
	if (val < 0) return ceil(val - 0.5);
	return floor(val + 0.5);
}

double step(Layer &layer)
{

	return -log(dis(gen)) / (layer.scattering + layer.absorption);
}
glm::dvec2 line_circle(double const x, double const y)
{
	
	auto const dx = x;
	auto const dy = y;
	auto const dr = sqrt(dx * dx + dy * dy);
	if(dr < radius)
	{
		return glm::dvec2(x, y);
	}
	
	auto const ix = (sign(dy)*dx*sqrt(radius*dr*dr)) / (dr*dr);
	auto const iy = (abs(dy)*sqrt(radius*dr*dr)) / (dr * dr);

	if(y < 0)
	{
		return glm::dvec2(-ix, -iy);

	}
	return glm::dvec2(ix, iy);

	
}
//Sampling recipe from the Improving Dwivedi Sampling Scheme paper
double dwivedi(Layer &layer)
{
	//TODO:
	//Problem: complex numbers are in the game
	auto const y = dis(gen);
	auto const x = (layer.v0 - 1) / (layer.v0 + 1);
	auto const wz_new = layer.v0*pow((layer.v0 + 1), y)*(pow((layer.v0 - 1), y) + pow((layer.v0 + 1), y));
	auto const tmp = pow(x, y);
	auto const wz = layer.v0 -(layer.v0 +1)*tmp;
	return (-log(1 - dis(gen))) / ((1 - wz / layer.v0)*(layer.scattering + layer.absorption));
}

void roulette(Photonstruct &photon)
{
	if (photon.weight <= 0.0)
	{
		photon.dead = true;
	}
	else if (dis(gen) < chance)
	{
		photon.weight *= chance;
	}
	else
	{
		photon.dead = true;
	}
}

void move(Photonstruct &photon, double const stepsize)
{
	
	photon.position.x = photon.position.x + photon.directions.x*stepsize;
	photon.position.y = photon.position.y + photon.directions.y*stepsize;
	photon.position.z = photon.position.z + photon.directions.z*stepsize;


	
	res.tissue.push_back(glm::vec3(photon.position.x, photon.position.y, photon.position.z));

}

glm::dvec2 calculate_scattering(material &material)
{
	double scattering;
	if (material.matproperties->anisotropy == 0.0)
	{	
		
		scattering = 2 * dis(gen) - 1;
	}
	else
	{
		scattering = 1 / 2 * material.matproperties->anisotropy
			* (1 + pow(material.matproperties->anisotropy, 2) - pow((1 - pow(material.matproperties->anisotropy, 2) /
																	(1 - material.matproperties->anisotropy +
																	2 * material.matproperties->anisotropy*dis(gen))), 2));
	}

	//Nachfragen ob gleiche random variable oder nicht
	auto const azimuthal = 2 * m_pi * dis(gen);

	return glm::dvec2(scattering, azimuthal);
} 

void update_direction(Photonstruct &photon, material &material)
{
	auto const angles = calculate_scattering(material);
	if (photon.directions.z > angle_threshold)
	{
		photon.directions.x = sin(acos(angles.x))*cos(angles.y);
		photon.directions.y = sin(acos(angles.x))*sin(angles.y);
		photon.directions.z = sin(sign(photon.directions.z))*angles.x;
	}
	else
	{
		photon.directions.x = sin(acos(angles.x)) / sqrt(1 - pow(photon.directions.z, 2))*(photon.directions.x*photon.directions.z*cos(angles.y) -
																				  photon.directions.y*sin(angles.y)) + photon.directions.x*angles.x;

		photon.directions.y = sin(acos(angles.x)) / sqrt(1 - pow(photon.directions.z, 2))*(photon.directions.y*photon.directions.z*cos(angles.y) +
																						  photon.directions.x*sin(angles.y)) + photon.directions.y*angles.x;

		photon.directions.z = -sin(acos(angles.x))*cos(angles.y)*sqrt(1 - pow(photon.directions.z, 2)) + photon.directions.z*angles.x;
	}
	
}

double specular_reflectance(Layer &layer)
{
	//First value is set to one because the specular reflectance is 1 of vacuum
	return pow((1 - layer.refrac), 2) / pow(1 + layer.refrac, 2);
}
double sign(double const x)
{
	return x > 0 ? 1 : -1;
}


double getv0(double const alpha)
{
	return alpha > 0.56 ? dvi.get_highv0(alpha) : dvi.get_lowv0(alpha);
}
bool is_hit(material &mat, Photonstruct &photon)
{
	//TODO: Test for both sampling sampling_scheme
	double stepsize;
	if(dwivedi_sam)
	{
		stepsize = dwivedi(*mat.matproperties);
	}
	else
	{
		stepsize = step(*mat.matproperties);
	}
	

	//TODO: Remove after tests
	//std::cout << stepsize << std::endl;
	photon.step = stepsize;
	if (photon.directions.z < 0)
	{
		//s1 = (z-z0) / (-mu_z) if mu_z < 0 :::::: z0 = 0
		auto const s1 = abs((photon.position.z) / photon.directions.z);
		//No check for lower boundery, because there is none
		if (s1 < stepsize)
		{
			move(photon, s1);
			stepsize -= s1;
			photon.step = stepsize;

			return true;
		}
	}

	return false;
}

void trace(material &mat, Photonstruct &photon)
{
	//Specular reflectance here, but not necassary because that would only decrement
	//the weight by constant


	if (is_hit(mat, photon))
	{
		const auto angleofincidence = acos(abs(photon.directions.z));
		const auto angleoftransmission = asin(mat.matproperties->refrac*sin(angleofincidence));


		const auto internalreflactance = 1 / 2 * ((pow(sin(angleofincidence - angleoftransmission), 2) / pow(sin(angleofincidence + angleoftransmission), 2) +
											  pow(tan(angleofincidence - angleoftransmission), 2) / pow(tan(angleofincidence + angleoftransmission), 2)));

		if (dis(gen) <= internalreflactance)
		{
			photon.directions.z *= -1;
		}
		else
		{
			res.transmittance.push_back(photon);
			updateBuckets(photon);
		}
	}
	else
	{
		
		move(photon, photon.step);
		absorption(mat, photon);
		update_direction(photon, mat);
		if(photon.weight < mat.wth && photon.weight > 0)
		{
			roulette(photon);
		}
		else if(photon.weight < 0)
		{
			photon.dead = true;
		}
		if (!photon.dead)
		{
			trace(mat, photon);
		}
	}
	
}

void render(material &mat)
{
	res.absorption = std::vector<glm::vec3>();
	res.transmittance = std::vector<Photonstruct>();
	res.tissue = std::vector<glm::vec3>();



	for (auto i = 0; i < mat.num_photons; i++)
	{
		Photonstruct photon;
		photon.dead = false;
		photon.directions = glm::vec3(0, 0, 1);
		photon.position = glm::vec3(0, 0, 0);
		photon.weight = 1;
		trace(mat, photon);
	}
}