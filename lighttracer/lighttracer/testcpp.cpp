
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

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

double dwivedi(double const scattering, double const absorption)
{
	//TODO:
	//Problem: complex numbers are in the game
	auto const y = dis(gen);
	auto const x = (dvi.highv0 - 1) / (dvi.highv0 + 1);


	auto const tmp = pow(x, y);
	auto const wz = dvi.highv0 - (dvi.highv0 + 1)*tmp;
	return (-log(1 - dis(gen))) / ((1 - wz / dvi.highv0)*(scattering + absorption));
}

int main()
{

	auto const absorption = 0.25;
	//Skin Scattering coefficient according to https://www.sciencedirect.com/science/article/pii/S0022202X15414836
	auto const scattering = 0.78;
	for(auto i = 0; i < 50; i++)
	{
		auto const res = dwivedi(scattering, absorption);
		std::cout << res << std::endl;
	}
	


}